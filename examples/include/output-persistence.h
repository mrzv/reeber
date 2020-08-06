#pragma once

#include <string>

#include <diy/master.hpp>
#include <diy/io/shared.hpp>
#include <dlog/log.h>

#include <reeber/triplet-merge-tree.h>

template<class Block, class LocalFunctor>
struct OutputPairs
{
    using RealType = typename Block::RealType;

    struct ExtraInfo
    {
        ExtraInfo(const std::string& _outfn, bool _verbose, diy::mpi::communicator& _world):
            outfn(_outfn), verbose(_verbose), world(_world),  ofs(_outfn, _world) {}

        std::string         outfn;
        bool                verbose;
        diy::mpi::communicator& world;
        diy::io::SharedOutFile ofs;
    };

    using Neighbor = typename Block::Neighbor;

    OutputPairs(const Block& _block, ExtraInfo* _extra, const LocalFunctor& _test_local, RealType _threshold, bool _ignore_zero_persistence):
        negate(_block.get_merge_tree().negate()),
        ignore_zero_persistence(_ignore_zero_persistence),
        threshold(_threshold),
        block(_block),
        extra(_extra),
        test_local(_test_local)
    {
    }

    void    operator()(Neighbor from, Neighbor through, Neighbor to) const
    {
        if (!test_local(block, from))
            return;

        if (extra->verbose)
        {
            std::string s;
            if (from != to)
                s = fmt::format("{} {} {} {} {} {}\n", from->vertex, from->value, through->vertex, through->value, to->vertex, to->value);
            else
                s = fmt::format("{} {} {} --\n",       from->vertex,  from->value, (block.get_merge_tree().negate() ? "-inf" : "inf"));
            extra->ofs << s;
        } else
        {
            RealType birth_time, death_time;
            birth_time = from->value;
            death_time = through->value;
            if (negate)
            {
//                LOG_SEV(info) << "in OutputPairs(), from = " << from->vertex << " recognized as local" << ", b = " << birth_time << ", d = " << death_time << ", to value = " << to->value << ", threshold = " << threshold;
                if (birth_time < threshold)
                    return;
                if (death_time < threshold)
                    //death_time = -std::numeric_limits<RealType>::infinity();
                    death_time = threshold;
            } else
            {
                if (birth_time > threshold)
                    return;
                if (death_time > threshold)
                    //death_time = std::numeric_limits<RealType>::infinity();
                    death_time = threshold;
            }

            if (ignore_zero_persistence and birth_time == death_time)
                return;

//            fmt::print("PERSISTENCE {} {}\n", birth_time, death_time);

            extra->ofs <<  birth_time << " " <<  death_time << "\n";
        }
    }

    const bool             negate;
    const bool             ignore_zero_persistence;
    const RealType         threshold;
    const Block&           block;
    ExtraInfo*             extra;
    const LocalFunctor&    test_local;
};

template<class Block, class LocalFunctor>
void output_persistence(Block* b, const diy::Master::ProxyWithLink& cp,
                        typename OutputPairs<Block, LocalFunctor>::ExtraInfo* extra,
                        const LocalFunctor& test_local,
                        typename Block::RealType threshold,
                        bool  _ignore_zero_persistence)
{
    if (b->get_merge_tree().size())
    {
//        LOG_SEV(info) << "Output persistence, block:   " << cp.gid() << ", vertices: " << b->get_merge_tree().n_vertices_total();
        reeber::traverse_persistence(b->get_merge_tree(),
                OutputPairs<Block, LocalFunctor>(*b, extra, test_local, threshold, _ignore_zero_persistence));
    }
}
