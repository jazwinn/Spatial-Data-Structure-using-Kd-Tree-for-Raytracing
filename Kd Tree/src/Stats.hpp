#ifndef STATS_HPP
#define STATS_HPP
#include <cstdint>
#include <vector>
namespace CS350 {
    struct Stats {
        unsigned            query_count              = 0;
        uint64_t            ray_vs_triangle          = 0;
        uint64_t            ray_vs_triangle_positive = 0;
        uint64_t            ray_vs_aabb              = 0;
        uint64_t            ray_vs_aabb_positive     = 0;
        uint64_t            kdtree_build_duration_ms = 0;
        std::vector<size_t> kdtree_last_traversed_nodes;

        void reset() {
            *this = {};
        }

        static Stats& instance() {
            static Stats instance;
            return instance;
        }
    };
}
#endif // STATS_HPP
