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
        ExtraInfo(const std::string& outfn_, bool verbose_, diy::mpi::communicator& _world):
            outfn(outfn_), verbose(verbose_), world(_world)   {}

        std::string         outfn;
        bool                verbose;
        diy::mpi::communicator& world;
    };

    using Neighbor = typename Block::Neighbor;

    OutputPairs(const Block& _block, const ExtraInfo& _extra, const LocalFunctor& _test_local, RealType _threshold, bool _ignore_zero_persistence):
        negate(_block.get_merge_tree().negate()),
        ignore_zero_persistence(_ignore_zero_persistence),
        threshold(_threshold),
        block(_block),
        extra(_extra),
        test_local(_test_local),
        ofs(extra.outfn, extra.world)
    {
    }

    void    operator()(Neighbor from, Neighbor through, Neighbor to) const
    {
        if (!test_local(block, from))
            return;

        if (extra.verbose)
        {
            std::string s;
            if (from != to)
                s = fmt::format("{} {} {} {} {} {}\n", from->vertex, from->value, through->vertex, through->value, to->vertex, to->value);
            else
                s = fmt::format("{} {} {} --\n",       from->vertex,  from->value, (block.get_merge_tree().negate() ? "-inf" : "inf"));
            ofs << s;
        } else
        {
            RealType birth_time, death_time;
            birth_time = from->value;
            death_time = through->value;
            if (negate)
            {
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
            ofs <<  birth_time <<  death_time << "\n";
        }
    }

    const bool             negate;
    const bool             ignore_zero_persistence;
    const RealType         threshold;
    const Block&           block;
    const ExtraInfo&       extra;
    const LocalFunctor&    test_local;
    mutable diy::io::SharedOutFile ofs;

};

template<class Block, class LocalFunctor>
void output_persistence(Block* b, const diy::Master::ProxyWithLink& cp,
                        const typename OutputPairs<Block, LocalFunctor>::ExtraInfo& extra,
                        const LocalFunctor& test_local,
                        typename Block::RealType threshold,
                        bool  _ignore_zero_persistence)
{
    LOG_SEV(debug) << "Block:   " << cp.gid();
    reeber::traverse_persistence(b->get_merge_tree(), OutputPairs<Block, LocalFunctor>(*b, extra, test_local, threshold, _ignore_zero_persistence));
}
