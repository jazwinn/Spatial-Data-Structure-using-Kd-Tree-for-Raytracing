#include "CS350Loader.hpp"
#include "Common.hpp"
#include "KdTree.hpp"
#include "Shapes.hpp"
#include "Stats.hpp"
#include "Utils.hpp"
#include "PRNG.h"

#include <fstream>
#include <limits>
#include <algorithm>

namespace {
    // Assets
    CS350::CS350PrimitiveData c_asset_bunny;
    CS350::CS350PrimitiveData c_asset_bunny_dense;
    CS350::CS350PrimitiveData c_asset_avocado;
    CS350::CS350PrimitiveData c_asset_dragon;

    struct Timer {
        using Clock = std::chrono::high_resolution_clock;
        Clock::time_point begin;
        Clock::time_point end;

        void start() {
            begin = std::chrono::high_resolution_clock::now();
        }

        auto stop() {
            end = std::chrono::high_resolution_clock::now();
            return duration_ms();
        }

        float duration_ms() const {
            auto duration = std::chrono::duration<float>(end - begin);
            return duration.count() * 1000.0f;
        }
    };
}

struct TestTriangle {
    CS350::Triangle geometry;
    size_t          index;

    operator CS350::Triangle() const { return geometry; };
    operator size_t() const { return index; };

    auto&       operator[](int i) { return geometry[i]; };
    auto const& operator[](int i) const { return geometry[i]; };
};

using Triangle = TestTriangle;
using Kdtree   = CS350::KdTree<TestTriangle>;

class Raytracing : public testing::Test {
  protected:

    static void SetUpTestSuite() {
        CS350::ChangeWorkdir();
        c_asset_bunny       = CS350::LoadCS350Binary("./assets/cs350/bunny.cs350_binary");
        c_asset_bunny_dense = CS350::LoadCS350Binary("./assets/cs350/bunny-dense.cs350_binary");
        c_asset_avocado     = CS350::LoadCS350Binary("./assets/cs350/avocado.cs350_binary");
        c_asset_dragon      = CS350::LoadCS350Binary("./assets/cs350/dragon.cs350_binary");
    }

    void SetUp() override {
        CS350::ChangeWorkdir();
        CS350::Stats::instance().reset();
    }

    static void add_mesh(std::vector<Triangle>&           triangles,
                         mat4 const&                      m2w,
                         CS350::CS350PrimitiveData const& primitive) {
        for (size_t i = 0; i < primitive.positions.size(); i += 3) {
            CS350::Triangle tri{};
            tri[0] = vec3(m2w * vec4(primitive.positions[i], 1));
            tri[1] = vec3(m2w * vec4(primitive.positions[i + 1], 1));
            tri[2] = vec3(m2w * vec4(primitive.positions[i + 2], 1));
            triangles.push_back({ tri, triangles.size() });
        }
    }

    static CS350::KdTreeConfig kdtreeconfig_basic() {
        CS350::KdTreeConfig cfg;
        cfg.cost_traversal    = 1.0f;
        cfg.cost_intersection = 80.0f;
        cfg.max_depth         = 20;
        cfg.min_triangles     = 20;
        return cfg;
    }

    static void save_graphs(Kdtree const& kdtree) {
        std::ofstream fdump("." + std::string(TestName()) + ".txt");
        kdtree.dump(fdump);
        std::string   graph_filename = "." + std::string(TestName()) + "-graph.dot";
        std::ofstream fgraph(graph_filename);
        kdtree.dump_graph(fgraph);
    }

    static bool is_lhs_contained_in_rhs(std::vector<Triangle> const& lhs,
                                        std::vector<Triangle> const& rhs) {
        for (auto lhs_element : lhs) {
            bool found_in_right = false;
            for (auto rhs_element : rhs) {
                if (lhs_element == rhs_element) {
                    found_in_right = true;
                    break;
                }
            }
            if (!found_in_right) {
                return true;
            }
        }
        return false;
    }

    static void ensure_node_sanity_internal(Kdtree const&            kdTree,
                                            size_t                   i,
                                            CS350::KdTreeNode const& node) {
        // Retrieve parent triangles (ALL of them)
        std::vector<Triangle> parent_tris;
        kdTree.get_triangles(i, parent_tris);

        // Retrieve triangles of left child (ALL of them)
        std::vector<Triangle> triangles_left;
        kdTree.get_triangles(i + 1, triangles_left);
        std::sort(triangles_left.begin(), triangles_left.end());
        // Retrieve triangles of right child (ALL of them)
        std::vector<Triangle> triangles_right;
        kdTree.get_triangles(node.next_child(), triangles_right);
        std::sort(triangles_right.begin(), triangles_right.end());
        ASSERT_LT(triangles_left.size(), parent_tris.size()) << "Child nodes should have less triangles than parent\n"
                                                             << "Problematic node: " << i;
        ASSERT_LT(triangles_right.size(), parent_tris.size()) << "Child nodes should have less triangles than parent\n"
                                                              << "Problematic node: " << i;
        ASSERT_FALSE(triangles_left.empty()) << "A child of an internal node should eventually contain triangles";
        ASSERT_FALSE(triangles_right.empty()) << "A child of an internal node should eventually contain triangles";
        ASSERT_TRUE(is_lhs_contained_in_rhs(triangles_right, triangles_left)) << "Node " << i + 1
                                                                              << " and " << node.next_child()
                                                                              << " are not disjoint, they have the same triangles";
        ASSERT_TRUE(is_lhs_contained_in_rhs(triangles_left, triangles_right)) << "Node " << i + 1
                                                                              << " and " << node.next_child()
                                                                              << " are not disjoint, they have the same triangles";
    }

    static void ensure_node_sanity_leaf(Kdtree const&            kdTree,
                                        size_t                   i,
                                        CS350::KdTreeNode const& node) {
        std::vector<Triangle> tris;
        kdTree.get_triangles(i, tris);
        ASSERT_EQ(node.primitive_count(), tris.size());
    }

    static void ensure_node_sanity(Kdtree const& kdTree, size_t i) {
        auto const& node = kdTree.nodes().at(i);

        if (node.is_internal()) {
            ASSERT_FALSE(node.is_leaf());
            ensure_node_sanity_internal(kdTree, i, node);
        } else {
            ASSERT_FALSE(node.is_internal());
            ASSERT_TRUE(node.is_leaf());
            ensure_node_sanity_leaf(kdTree, i, node);
        }
    }

    static void ensure_nodes_sanity(Kdtree const& kdTree) {
        for (size_t i = 0; i < kdTree.nodes().size(); ++i) {
            ensure_node_sanity(kdTree, i);
        }
    }

    static void TestGenericBuild(CS350::CS350PrimitiveData const& asset,
                                 int                              max_depth = 20) {
        Kdtree kdtree;
        kdtree.cfg()           = kdtreeconfig_basic();
        kdtree.cfg().max_depth = max_depth;

        std::vector<Triangle> triangles;
        add_mesh(triangles, mat4(1), asset);

        Timer timer_build;
        timer_build.start();
        kdtree.build(triangles);
        std::cout << fmt::format("Build took: {:02f}ms", timer_build.stop()) << std::endl;
        std::cout << fmt::format("Kdtree depth: {}", kdtree.get_depth()) << std::endl;
        std::cout << fmt::format("Kdtree nodes: {}", kdtree.nodes().size()) << std::endl;
        std::cout << fmt::format("Kdtree triangles: {}", kdtree.triangles().size()) << std::endl;

        save_graphs(kdtree);
        ensure_nodes_sanity(kdtree);
    }

    static CS350::Ray::Intersection raytest_bruteforce(std::vector<Triangle> const& triangles,
                                                       CS350::Ray const&            ray) {
        CS350::Ray::Intersection best_intersection = { -1.0f };
        for (auto const& tri : triangles) {
            CS350::Stats::instance().ray_vs_triangle++;
            auto intersection = ray.intersect(tri);
            if (intersection) {
                // Debug
                // std::cerr << fmt::format("[BruteForce] Ray {{{}}} intersected {{{}}}", ray, tri.geometry) << std::endl;
                CS350::Stats::instance().ray_vs_triangle_positive++;
                if (intersection.t < best_intersection.t || (!best_intersection)) {
                    best_intersection = intersection.t;
                }
            }
        }
        return best_intersection;
    }

    static CS350::Ray::Intersection raytest_kdtree(Kdtree const&     kdtree,
                                                   CS350::Ray const& ray) {
        CS350::Ray::Intersection best_intersection = { -1.0f };
        kdtree.get_closest(ray, best_intersection.t);
        return best_intersection;
    }

    static void test_generic_ray_test(CS350::CS350PrimitiveData const& asset,
                                      mat4 const&                      m2w        = mat4(1),
                                      int                              max_depth  = 20,
                                      int                              test_count = 10000) {
        Kdtree kdtree;
        kdtree.cfg()           = kdtreeconfig_basic();
        kdtree.cfg().max_depth = max_depth;

        std::vector<Triangle> triangles;
        add_mesh(triangles, m2w, asset);
        kdtree.build(triangles);
        save_graphs(kdtree);

        CS350::Stats averages_brute_force{};
        CS350::Stats averages_kdtree{};
        for (int i = 0; i < test_count; i++) {
            // Create a random ray
            CS350::Ray ray{};
            ray.start.x = CS170::Utils::Random(-50.0f, 50.0f);
            ray.start.y = CS170::Utils::Random(-50.0f, 50.0f);
            ray.start.z = CS170::Utils::Random(-50.0f, 50.0f);
            vec3 target;
            target.x = CS170::Utils::Random(-1.0f, 1.0f);
            target.y = CS170::Utils::Random(-1.0f, 1.0f);
            target.z = CS170::Utils::Random(-1.0f, 1.0f);
            target   = {};
            ray.dir  = target - ray.start;

            // Brute force
            CS350::Stats::instance().reset();
            auto bruteforce_solution = raytest_bruteforce(triangles, ray);
            auto bruteforce_stats    = CS350::Stats::instance();

            // KdTree
            CS350::Stats::instance().reset();
            auto kdtree_solution = raytest_kdtree(kdtree, ray);
            auto kdtree_stats    = CS350::Stats::instance();

            // Stats keep track
            averages_brute_force.ray_vs_triangle += bruteforce_stats.ray_vs_triangle;
            averages_brute_force.ray_vs_triangle_positive += bruteforce_stats.ray_vs_triangle_positive;
            averages_brute_force.ray_vs_aabb += bruteforce_stats.ray_vs_aabb;
            averages_brute_force.ray_vs_aabb_positive += bruteforce_stats.ray_vs_aabb_positive;

            averages_kdtree.ray_vs_triangle += kdtree_stats.ray_vs_triangle;
            averages_kdtree.ray_vs_triangle_positive += kdtree_stats.ray_vs_triangle_positive;
            averages_kdtree.ray_vs_aabb += kdtree_stats.ray_vs_aabb;
            averages_kdtree.ray_vs_aabb_positive += kdtree_stats.ray_vs_aabb_positive;

            // Debug
            // std::cout << "=======\n";
            // std::cout << fmt::format("kdtree_stats.kdtree_build_duration_ms: {}", kdtree_stats.kdtree_build_duration_ms) << std::endl;
            // std::cout << fmt::format("kdtree_stats.query_count: {}", kdtree_stats.query_count) << std::endl;
            // std::cout << fmt::format("kdtree_stats.ray_vs_aabb: {}", kdtree_stats.ray_vs_aabb) << std::endl;
            // std::cout << fmt::format("kdtree_stats.ray_vs_aabb_positive: {}", kdtree_stats.ray_vs_aabb_positive) << std::endl;
            // std::cout << fmt::format("kdtree_stats.ray_vs_triangle: {}", kdtree_stats.ray_vs_triangle) << std::endl;
            // std::cout << fmt::format("kdtree_stats.ray_vs_triangle_positive: {}", kdtree_stats.ray_vs_triangle_positive) << std::endl;

            // Ensure same results (with KdTree being faster)
            ASSERT_FLOAT_EQ(kdtree_solution.t, bruteforce_solution.t)
                << "Bruteforce solution and KdTree solution differ!" << "\n"
                << "with ray: " << ray << "\n"
                << "iteration: " << i;
            ASSERT_LT(kdtree_stats.ray_vs_triangle, bruteforce_stats.ray_vs_triangle)
                << "KdTree is not faster than brute force" << "\n"
                << "with ray: " << ray << "\n"
                << "iteration: " << i;

            // Check traversed nodes
            auto const& traversed_nodes = CS350::Stats::instance().kdtree_last_traversed_nodes;
            if (kdtree_solution) {
                ASSERT_FALSE(traversed_nodes.empty()) << "Traversed nodes must be recorded into the stats system";
            }

            // Must have intersected
            for (auto const& node : traversed_nodes) {
                auto const& aabb             = kdtree.aabbs().at(node);
                auto        ray_intersection = ray.intersect(aabb);

                // Intersection must exist
                ASSERT_TRUE(ray_intersection) << "A node is traversed that is not intersected by the ray." << "\n"
                                              << "with ray: " << ray << "\n"
                                              << "with aabb: " << aabb << "\n"
                                              << "iteration: " << i;
            }

            // Ensure that the KdTree is 'good'
            ASSERT_LT(kdtree_stats.ray_vs_triangle, bruteforce_stats.ray_vs_triangle / 6)
                << "KdTree is not fast enough" << "\n"
                << "with ray: " << ray << "\n"
                << "iteration: " << i;
        }

        // Averages check

        float averages_brute_force_ray_vs_triangle = (float)(averages_brute_force.ray_vs_triangle) / (float)test_count;
        float averages_kdtree_ray_vs_triangle      = (float)(averages_kdtree.ray_vs_triangle) / (float)test_count;
        float averages_kdtree_ray_vs_aabb          = (float)(averages_kdtree.ray_vs_aabb) / (float)test_count;

        // Debug
        std::cout << fmt::format("averages_brute_force_ray_vs_triangle: {}", averages_brute_force_ray_vs_triangle) << std::endl;
        std::cout << fmt::format("averages_kdtree_ray_vs_triangle: {}", averages_kdtree_ray_vs_triangle) << std::endl;
        std::cout << fmt::format("averages_kdtree_ray_vs_aabb: {}", averages_kdtree_ray_vs_aabb) << std::endl;
    }
};

//// Build only
TEST_F(Raytracing, Build_Avocado) { TestGenericBuild(c_asset_avocado); }
TEST_F(Raytracing, Build_Bunny8) { TestGenericBuild(c_asset_bunny, 8); }
TEST_F(Raytracing, Build_Bunny16) { TestGenericBuild(c_asset_bunny, 16); }
TEST_F(Raytracing, Build_BunnyDense8) { TestGenericBuild(c_asset_bunny_dense, 8); }
TEST_F(Raytracing, Build_BunnyDense16) { TestGenericBuild(c_asset_bunny_dense, 16); }
TEST_F(Raytracing, Build_BunnyDense32) { TestGenericBuild(c_asset_bunny_dense, 32); }
TEST_F(Raytracing, Build_Dragon8) { TestGenericBuild(c_asset_dragon, 8); }
TEST_F(Raytracing, Build_Dragon16) { TestGenericBuild(c_asset_dragon, 16); }
TEST_F(Raytracing, Build_Dragon32) { TestGenericBuild(c_asset_dragon, 32); }
//
//// Raytesting
TEST_F(Raytracing, Raytest_Bunny8) { test_generic_ray_test(c_asset_bunny, glm::translate(vec3(0, 0, 2)), 8); }
TEST_F(Raytracing, Raytest_Bunny16) { test_generic_ray_test(c_asset_bunny, glm::translate(vec3(0, 0, 2)), 16); }
TEST_F(Raytracing, Raytest_BunnyDense8) { test_generic_ray_test(c_asset_bunny_dense, glm::translate(vec3(0, 0, 2)), 8); }
TEST_F(Raytracing, Raytest_BunnyDense16) { test_generic_ray_test(c_asset_bunny_dense, glm::translate(vec3(0, 0, 2)), 16); }
TEST_F(Raytracing, Raytest_BunnyDense32) { test_generic_ray_test(c_asset_bunny_dense, glm::translate(vec3(0, 0, 2)), 32); }
TEST_F(Raytracing, Raytest_Dragon8) { test_generic_ray_test(c_asset_dragon, mat4(1), 8, 1000); }
TEST_F(Raytracing, Raytest_Dragon16) { test_generic_ray_test(c_asset_dragon, mat4(1), 16, 1000); }
TEST_F(Raytracing, Raytest_Dragon32) { test_generic_ray_test(c_asset_dragon, mat4(1), 32, 1000); }

 //Manual tests
TEST_F(Raytracing, ManualNormal) {
    // Every object is going to have an AABB that is 1x1 (consider it clsoe to a 2D problem)
    std::vector<vec3> objects_bottom_left_corner;
    objects_bottom_left_corner.push_back({});
    objects_bottom_left_corner.push_back({ 1, 2, 0 });
    objects_bottom_left_corner.push_back({ 2, 3, 0 });
    objects_bottom_left_corner.push_back({ 4, 0, 0 });
    objects_bottom_left_corner.push_back({ 6, 0, 0 });
    objects_bottom_left_corner.push_back({ 7, 2, 0 });

    // Convert it to triangles so it's processed by the kdtree
    std::vector<Triangle> triangles;
    triangles.reserve(objects_bottom_left_corner.size());
    for (auto& corner : objects_bottom_left_corner) {
        triangles.push_back({ CS350::Triangle(corner, corner + vec3(1, 0, 0), corner + vec3(1, 1, 0)), 0 });
    }

    // Create KdTree
    Kdtree kdtree;
    kdtree.cfg().cost_intersection = 80;   //
    kdtree.cfg().cost_traversal    = 1;    //
    kdtree.cfg().max_depth         = 1000; // Irrelevant
    kdtree.cfg().min_triangles     = 3;    // Don't split if fewer than min
    kdtree.build(triangles);
    save_graphs(kdtree);

    //
    ASSERT_EQ(kdtree.get_depth(), 1) << "Node should be split once";
    ASSERT_EQ(kdtree.nodes().size(), 3) << "Incorrect amount of nodes";
    ASSERT_TRUE(kdtree.nodes()[0].is_internal()) << "Root should be internal";
    ASSERT_TRUE(kdtree.nodes()[1].is_leaf()) << "Other nodes should be leaf";
    ASSERT_EQ(kdtree.nodes()[1].primitive_count(), 3) << "Bad split";
    ASSERT_TRUE(kdtree.nodes()[2].is_leaf()) << "Other nodes should be leaf";
    ASSERT_TRUE(kdtree.nodes()[2].is_leaf()) << "Other nodes should be leaf";
    ASSERT_EQ(kdtree.nodes()[2].primitive_count(), 3) << "Bad split";
}

TEST_F(Raytracing, ManualLeaf) {
    // Every object is going to have an AABB that is 1x1 (consider it clsoe to a 2D problem)
    std::vector<vec3> objects_bottom_left_corner;
    objects_bottom_left_corner.push_back({});
    objects_bottom_left_corner.push_back({ 1, 2, 0 });
    objects_bottom_left_corner.push_back({ 2, 3, 0 });
    objects_bottom_left_corner.push_back({ 4, 0, 0 });
    objects_bottom_left_corner.push_back({ 6, 0, 0 });
    objects_bottom_left_corner.push_back({ 7, 2, 0 });

    // Convert it to triangles so it's processed by the kdtree
    std::vector<Triangle> triangles;
    triangles.reserve(objects_bottom_left_corner.size());
    for (auto& corner : objects_bottom_left_corner) {
        triangles.push_back({ CS350::Triangle(corner, corner + vec3(1, 0, 0), corner + vec3(1, 1, 0)), 0 });
    }

    // Create KdTree
    Kdtree kdtree;
    kdtree.cfg().cost_intersection = 64;   // Keep it "low"
    kdtree.cfg().cost_traversal    = 1000; // Making this a high cost, favours making a leaf
    kdtree.cfg().max_depth         = 1000; // Irrelevant
    kdtree.cfg().min_triangles     = 1;    // Don't delay creation
    kdtree.build(triangles);
    save_graphs(kdtree);

    //
    ASSERT_EQ(kdtree.get_depth(), 0) << "Best option is to make a leaf!";
    ASSERT_TRUE(kdtree.nodes().front().is_leaf());
    ASSERT_EQ(kdtree.nodes().front().primitive_count(), triangles.size());
}
