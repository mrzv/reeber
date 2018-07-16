#pragma once

#include <string>

#include <diy/master.hpp>
#include <dlog/log.h>

#include <reeber/triplet-merge-tree.h>
#include <reeber/triplet-merge-tree.h>

template<class Block, class LocalFunctor>
struct OutputPairs
{
    struct ExtraInfo
    {
        ExtraInfo(const std::string& outfn_, bool verbose_):
            outfn(outfn_), verbose(verbose_)   {}

        std::string         outfn;
        bool                verbose;
    };

    using Neighbor = typename Block::TripletMergeTree::Neighbor;

    OutputPairs(const Block& block_, const ExtraInfo& extra_, const LocalFunctor& test_local_):
        block(block_),
        extra(extra_),
        test_local(test_local_)
    {
        std::string   dgm_fn = fmt::format("{}-b{}.dgm", extra.outfn, block.gid);
        ofs.open(dgm_fn.c_str());
    }

    void    operator()(Neighbor from, Neighbor through, Neighbor to) const
    {
        if (!test_local(block, from))
            return;

        if (extra.verbose)
        {
            if (from != to)
                fmt::print(ofs, "{} {} {} {} {} {}\n", from->vertex, from->value, through->vertex, through->value, to->vertex, to->value);
            else
                fmt::print(ofs, "{} {} {} --\n",       from->vertex,  from->value, (block.get_merge_tree().negate() ? "-inf" : "inf"));
        } else
            fmt::print(ofs, "{} {}\n", from->value, through->value);
    }

    const Block&           block;
    const ExtraInfo&       extra;
    const LocalFunctor&    test_local;
    mutable std::ofstream  ofs;
};

template<class Block, class LocalFunctor>
void output_persistence(Block* b, const diy::Master::ProxyWithLink& cp,
                        const typename OutputPairs<Block, LocalFunctor>::ExtraInfo& extra,
                        const LocalFunctor& test_local)
{
    LOG_SEV(debug) << "Block:   " << cp.gid();
    reeber::traverse_persistence(b->get_merge_tree(), OutputPairs<Block, LocalFunctor>(*b, extra, test_local));
}
