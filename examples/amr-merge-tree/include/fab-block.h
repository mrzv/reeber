#pragma once

#include <diy/serialization.hpp>
#include <diy/grid.hpp>
#include <diy/vertices.hpp>
#include <diy/fmt/format.h>
#include <diy/fmt/ostream.h>


template<class T, unsigned D>
struct FabBlock
{
    using Shape = diy::Point<int, D>;
    using Vertex = diy::Point<int, D>;
    using Grid = diy::Grid<T, D>;
    using GridRef = diy::GridRef<T, D>;

    FabBlock() :
            fab(fab_storage_.data(), fab_storage_.shape(), fab_storage_.c_order())
    {
    }

    FabBlock(T* data, const Shape& shape) :
            fab(data, shape, /* c_order = */ false)
    {
    }

    static void* create()
    {
        return new FabBlock;
    }

    static void destroy(void* b)
    {
        delete static_cast<FabBlock*>(b);
    }

    static void save(const void* b_, diy::BinaryBuffer& bb);

    static void load(void* b_, diy::BinaryBuffer& bb);

    diy::Grid<T, D> fab_storage_;        // container, in case we own the data
    diy::GridRef<T, D> fab;

};

template<class T, unsigned D>
void
FabBlock<T, D>::save(const void* b_, diy::BinaryBuffer& bb)
{
    auto* b = static_cast<const FabBlock<T, D>*>(b_);
    diy::save(bb, b->fab.shape());
    diy::save(bb, b->fab.c_order());
    diy::save(bb, b->fab.data(), b->fab.size());
}

template<class T, unsigned D>
void
FabBlock<T, D>::load(void* b_, diy::BinaryBuffer& bb)
{
    auto* b = static_cast<FabBlock<T, D>*>(b_);

    Shape shape;
    bool c_order;
    diy::load(bb, shape);
    diy::load(bb, c_order);

    b->fab_storage_ = decltype(b->fab_storage_)(shape, c_order);
    diy::load(bb, b->fab_storage_.data(), b->fab_storage_.size());

    b->fab = decltype(b->fab)(b->fab_storage_.data(), shape, c_order);     // fab points to the data in fab_storage_
}

template<class T, unsigned D>
void change_to_c_order(FabBlock<T, D>* b)
{
    using Block = FabBlock<T, D>;
    using Vertex = typename Block::Vertex;
    using Grid = typename Block::Grid;

    if (b->fab_storage_.c_order())
        return;

    Grid tmp_grid(b->fab_storage_.shape(), true);

    diy::for_each(b->fab_storage_.shape(), [&](const Vertex v) {
        tmp_grid(v) = b->fab(v);
    });

    b->fab_storage_.swap(tmp_grid);

    assert(b->fab_storage_.c_order());

    b->fab = decltype(b->fab)(b->fab_storage_.data(), b->fab_storage_.shape(),
                              b->fab_storage_.c_order());     // fab points to the data in fab_storage_

}
