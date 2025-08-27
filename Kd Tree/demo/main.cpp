#include "DebugRenderer.hpp"
#include "Primitive.hpp"
#include "Utils.hpp"
#include "OpenGl.hpp"
#include "Math.hpp"
#include "Window.hpp"
#include "Camera.hpp"
#include "ImGui.hpp"
#include "CS350Loader.hpp"
#include "Stats.hpp"

#include <PRNG.h>
#include <imgui.h>
#include <memory>
#include <queue>
#include <thread>

#include "../src/KdTree.hpp"

namespace {
    // Assets
    char const *const c_asset_bunny = "./assets/cs350/bunny.cs350_binary";
    char const *const c_asset_bunny_dense = "./assets/cs350/bunny-dense.cs350_binary";
    char const *const c_asset_avocado = "./assets/cs350/avocado.cs350_binary";
    char const *const c_asset_cube = "./assets/cs350/cube.cs350_binary";
    char const *const c_asset_dragon = "./assets/cs350/dragon.cs350_binary";

    // Camera controls
    float cCameraSpeed = 0.2f;

    void ScrollCallback(GLFWwindow * /*window*/, double /*xoffset*/, double yoffset) {
        if (!ImGui::GetIO().WantCaptureMouse) {
            cCameraSpeed += (float) yoffset * 0.05f;
            cCameraSpeed = glm::clamp(cCameraSpeed, 0.0f, 1.0f);
        }
    }

    void CameraMovementFly(CS350::Camera &camera, float dt, CS350::Window &window) {
        auto *glfwWindow = window.handle();
        auto camera_dir = camera.target() - camera.position();
        auto camera_position = camera.position();
        static vec2 cursor_position{};

        // Cursor
        double cursor_x = 0.0;
        double cursor_y = 0.0;
        glfwGetCursorPos(glfwWindow, &cursor_x, &cursor_y); {
            // Input
            if (glfwGetMouseButton(glfwWindow, GLFW_MOUSE_BUTTON_2) != 0) {
                float speed = glm::mix(0.01f, 50.0f, cCameraSpeed * cCameraSpeed);
                auto side = normalize(cross(camera_dir, {0, 1, 0}));

                if (glfwGetKey(glfwWindow, GLFW_KEY_LEFT_SHIFT)) {
                    speed *= 4.0f;
                }
                if (glfwGetKey(glfwWindow, GLFW_KEY_LEFT_ALT)) {
                    speed /= 4.0f;
                }
                // Move
                if (glfwGetKey(glfwWindow, GLFW_KEY_W)) {
                    camera_position += normalize(camera_dir) * dt * speed;
                }
                if (glfwGetKey(glfwWindow, GLFW_KEY_S)) {
                    camera_position -= normalize(camera_dir) * dt * speed;
                }
                if (glfwGetKey(glfwWindow, GLFW_KEY_A)) {
                    camera_position -= normalize(side) * dt * speed;
                }
                if (glfwGetKey(glfwWindow, GLFW_KEY_D)) {
                    camera_position += normalize(side) * dt * speed;
                }

                // View
                vec2 cursor_delta = {(float) cursor_x - cursor_position.x, (float) cursor_y - cursor_position.y};
                camera_dir = vec3(vec4(camera_dir, 0) * rotate(glm::radians(15.0f) * 0.01f * cursor_delta.y, side));
                camera_dir = vec3(
                    vec4(camera_dir, 0) * rotate(glm::radians(15.0f) * 0.01f * cursor_delta.x, vec3(0, 1, 0)));
            }
            cursor_position = {(float) cursor_x, (float) cursor_y};
            auto camera_target = camera_position + camera_dir;

            // Submit changes
            camera.SetPosition(camera_position);
            camera.SetTarget(camera_target);
            camera.SetProjection(camera.fov_deg(), window.size(), camera.near(), camera.far());
        }
    }
}

struct GameObject {
    std::shared_ptr<CS350::Primitive> primitive;
    CS350::CS350PrimitiveData primitive_data;

    // m2w
    vec3 position;
    vec3 scale;
    quat rotation;
    mat4 m2w;

    void compute_m2w() {
        m2w = glm::translate(position) * glm::toMat4(rotation) * glm::scale(scale);
    }
};

struct Scene {
    std::vector<GameObject> objects;

    // Raytracing specific
    CS350::KdTree<CS350::Triangle> kdtree;

    Scene() {
    }

    ~Scene() {
    }

    void add_object(char const *filename) {
        GameObject obj;

        // Load CS350
        obj.primitive_data = CS350::LoadCS350Binary(filename);
        std::vector<float> vbo;
        vbo.reserve(obj.primitive_data.positions.size() * 3);
        for (auto pos: obj.primitive_data.positions) {
            vbo.push_back(pos.x);
            vbo.push_back(pos.y);
            vbo.push_back(pos.z);
        }

        // Create OpenGL asset
        obj.primitive = std::make_shared<CS350::Primitive>(vbo.data(), vbo.size() * sizeof(float));

        // m2w
        obj.position = {};
        obj.scale = {1, 1, 1};
        obj.rotation = glm::normalize(glm::quat(vec3(glm::radians(00.0f), 0, 0)));
        obj.compute_m2w();
        objects.push_back(obj);
    }

    void reconstruct_scene() {
        // Prepare scene
        auto cfg = kdtree.cfg();
        std::vector<CS350::Triangle> triangles;
        for (auto const &obj: objects) {
            // Add geometry
            for (size_t i = 0; i < obj.primitive_data.positions.size(); i += 3) {
                CS350::Triangle tri{};
                tri[0] = vec3(obj.m2w * vec4(obj.primitive_data.positions[i], 1));
                tri[1] = vec3(obj.m2w * vec4(obj.primitive_data.positions[i + 1], 1));
                tri[2] = vec3(obj.m2w * vec4(obj.primitive_data.positions[i + 2], 1));
                triangles.push_back(tri);
            }
        }
        kdtree.build(triangles);
    }
};


struct {
    bool draw_triangles = false;
    float kdtree_aabb_alpha = 0.025f;
    CS350::Ray last_manual_ray = {vec3(std::numeric_limits<float>::quiet_NaN()), vec3()};
    bool debug_manual_ray = false;
    bool debug_kdtree = false;
    int debug_manual_ray_traversed_node = -1;
    std::vector<size_t> path;
    std::string node_info;

    vec4 node_color(int node_index) const {
        srand(node_index + 20);
        return vec4(glm::linearRand(0.1f, 1.0f), glm::linearRand(0.1f, 1.0f), glm::linearRand(0.1f, 1.0f),
                    kdtree_aabb_alpha);
    }
} debug_options;

int main() {
    CS350::ChangeWorkdir("bin");
    CS350::Window::initialize_system();
    CS350::Window w({1920, 1080});
    glfwSetScrollCallback(w.handle(), ScrollCallback);
    ImGuiInitialize(w.handle()); {
        //
        CS350::DebugRenderer debug;
        Scene scene;
        CS350::Camera camera;
        camera.SetPosition({0, 0, 2});
        camera.SetTarget({0, 0, 0});
        camera.set_fov_deg(50.0f); {
            // Sample scene

            // Object
            // scene.add_object(c_asset_avocado);
            // scene.add_object(c_asset_cube);
            scene.add_object(c_asset_bunny);

            // Test values
            auto &cfg = scene.kdtree.cfg();
            cfg.cost_traversal = 1.0f;
            cfg.cost_intersection = 80.0f;
            cfg.max_depth = 20;
            cfg.min_triangles = 20; {
                // Manual ray
                auto &ray = debug_options.last_manual_ray;
                ray.start = {-0.973711, -1.52317, -3.02257};
                ray.dir = {-0.217702, 0.329, -0.127437};
            }
            scene.reconstruct_scene();
        }

        //
        while (!w.should_exit()) {
            // CS170::Utils::srand(5, 5);
            float dt = ImGui::GetIO().DeltaTime;
            // Window
            w.update();

            // Camera matrices
            CameraMovementFly(camera, dt, w);
            camera.Update();
            mat4 v = camera.view();
            mat4 p = camera.proj();
            mat4 vp = camera.viewProj();

            // General rendering states
            ivec2 windowSize = w.size();
            glDepthMask(GL_TRUE); // Restore
            glViewport(0, 0, windowSize.x, windowSize.y);
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            glClearDepth(1.0);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LEQUAL);
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE);
            glDepthMask(GL_FALSE);

            glDisable(GL_CULL_FACE);
            debug.draw_begin();
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDisable(GL_DEPTH_TEST);
            glDepthMask(GL_FALSE);
            for (auto &obj: scene.objects) {
                debug.draw_primitive(vp, *obj.primitive, obj.m2w, vec4(0.25f, 0.25f, 0.25f, 0.15f));
            }
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

            // For debugging
            auto &kdtree = scene.kdtree;
            auto debug_draw_node = [&](int node_index) {
                auto const &n = kdtree.nodes().at(node_index);

                // Aabb
                glEnable(GL_BLEND);
                glDisable(GL_CULL_FACE);
                glDepthMask(GL_FALSE);
                glDisable(GL_DEPTH_TEST);
                auto color = debug_options.node_color(node_index);
                auto aabb = kdtree.aabbs().at(node_index);
                debug.draw_aabb(vp, aabb.get_center(), aabb.get_extents(), color);
                debug.draw_aabb_wireframe(vp, aabb.get_center(), aabb.get_extents(), color);
                // Split
                if (n.is_internal()) {
                    aabb.min[n.axis()] = aabb.max[n.axis()] = n.split();
                    debug.draw_aabb(vp, aabb.get_center(), aabb.get_extents(), color);
                    debug.draw_aabb_wireframe(vp, aabb.get_center(), aabb.get_extents(), color);
                }

                if (debug_options.draw_triangles) {
                    // Triangles
                    glDisable(GL_BLEND);
                    glDisable(GL_CULL_FACE);
                    glDepthMask(GL_TRUE);
                    glEnable(GL_DEPTH_TEST);
                    std::vector<CS350::Triangle> triangles;
                    kdtree.get_triangles(node_index, triangles);
                    for (auto &t: triangles) {
                        debug.draw_triangle(vp, t[0], t[1], t[2], color);
                        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                        debug.draw_triangle(vp, t[0], t[1], t[2], {0, 0, 0, 1});
                        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    }
                }
            };
            ImGuiNewFrame();

            if (ImGui::Begin("KdTree")) {
                if (ImGui::Button("Avocado")) {
                    scene.objects.clear();
                    scene.add_object(c_asset_avocado);
                    scene.reconstruct_scene();
                }
                if (ImGui::Button("Bunny")) {
                    scene.objects.clear();
                    scene.add_object(c_asset_bunny);
                    scene.reconstruct_scene();
                }
                if (ImGui::Button("Bunny dense")) {
                    scene.objects.clear();
                    scene.add_object(c_asset_bunny_dense);
                    scene.reconstruct_scene();
                }
                if (ImGui::Button("Cube")) {
                    scene.objects.clear();
                    scene.add_object(c_asset_cube);
                    scene.reconstruct_scene();
                }
                if (ImGui::Button("Dragon")) {
                    scene.objects.clear();
                    scene.add_object(c_asset_dragon);
                    scene.reconstruct_scene();
                }

                ImGui::Text("ray_vs_aabb: %lu", CS350::Stats::instance().ray_vs_aabb);
                ImGui::Text("ray_vs_aabb_positive: %lu", CS350::Stats::instance().ray_vs_aabb_positive);
                ImGui::Text("ray_vs_triangle: %lu", CS350::Stats::instance().ray_vs_triangle);
                ImGui::Text("ray_vs_triangle_positive: %lu", CS350::Stats::instance().ray_vs_triangle_positive);


                ImGui::DragFloat("cost_intersection", &kdtree.cfg().cost_intersection, 0.5f, 0.0f, FLT_MAX);
                ImGui::DragFloat("cost_traversal", &kdtree.cfg().cost_traversal, 0.5f, 0.0f, FLT_MAX);
                ImGui::DragInt("max_depth", &kdtree.cfg().max_depth, 0.1f, 0, INT_MAX);
                ImGui::DragInt("min_triangles", &kdtree.cfg().min_triangles, 0.1f, 1, INT_MAX);
                if (ImGui::Button("Reconstruct")) {
                    scene.reconstruct_scene();
                }

                ImGui::Text("Node count: %d", scene.kdtree.nodes().size());
                ImGui::Text("Triangle count: %d", scene.kdtree.triangles().size());
                ImGui::Text("Depth: %d", scene.kdtree.get_depth());
                ImGui::Separator();
                ImGui::Text("DEBUG");
                ImGui::SliderFloat("Debug alpha", &debug_options.kdtree_aabb_alpha, 0.0f, 1.0f);
                ImGui::Checkbox("Draw triangles recursively", &debug_options.draw_triangles);
                if (ImGui::Checkbox("Debug ray", &debug_options.debug_manual_ray)) {
                    if (debug_options.debug_manual_ray) {
                        debug_options.debug_kdtree = false;
                    }
                }
                if (ImGui::Checkbox("Debug nodes", &debug_options.debug_kdtree)) {
                    if (debug_options.debug_kdtree) {
                        debug_options.debug_manual_ray = false;
                    }
                }

                // Manual ray
                if (debug_options.debug_manual_ray) {
                    auto &ray = debug_options.last_manual_ray;
                    ImGui::DragFloat3("Debug ray start", &ray.start.x);
                    ImGui::DragFloat3("Debug ray dir", &ray.dir.x);
                    vec3 end = ray.start + ray.dir;
                    if (ImGui::Guizmo(&ray.start, v, p)) {
                        ray.dir = end - ray.start;
                    }
                    if (ImGui::Guizmo(&end, v, p)) {
                        ray.dir = end - ray.start;
                    }
                    ImGui::SliderInt("Debug node", &debug_options.debug_manual_ray_traversed_node, -1,
                                     CS350::Stats::instance().kdtree_last_traversed_nodes.size() - 1);
                }

                // Draw manual ray nodes
                if (debug_options.debug_manual_ray && debug_options.debug_manual_ray_traversed_node < int(
                        CS350::Stats::instance().kdtree_last_traversed_nodes.size())) {
                    auto const &nodes = CS350::Stats::instance().kdtree_last_traversed_nodes;
                    if (debug_options.debug_manual_ray_traversed_node == -1) {
                        for (auto &node_index: nodes) {
                            debug_draw_node(node_index);
                        }
                    } else {
                        auto node_index = nodes[debug_options.debug_manual_ray_traversed_node];
                        debug_draw_node(node_index);

                        // Debug intersections
                        auto const &node = kdtree.nodes()[node_index];
                        if (node.is_internal()) {
                            auto axis = node.axis();
                            auto const &ray = debug_options.last_manual_ray;
                            auto const &bv = kdtree.aabbs().at(node_index);
                            float entry_t = 0.0f;
                            float exit_t = FLT_MAX;
                            for (int axis = 0; axis < 3; ++axis) {
                                auto entry = (bv.min[axis] - ray.start[axis]) / ray.dir[axis];
                                auto exit = (bv.max[axis] - ray.start[axis]) / ray.dir[axis];
                                if (ray.dir[axis] < 0) {
                                    std::swap(entry, exit);
                                }
                                entry_t = glm::max(entry_t, entry);
                                exit_t = glm::min(exit_t, exit);
                            }
                            auto split = node.split();
                            auto node_t = (split - ray.start[axis]) / ray.dir[axis];
                            debug.draw_sphere(vp, camera.position(), ray.evaluate(entry_t), 0.05f,
                                              {0.25, 0.9, 0.1, 1.0});
                            debug.draw_sphere(vp, camera.position(), ray.evaluate(exit_t), 0.05f,
                                              {0.9, 0.25, 0.1, 1.0});
                            debug.draw_sphere(vp, camera.position(), ray.evaluate(node_t), 0.05f,
                                              {0.9, 0.25, 0.9, 1.0});
                        }
                    }
                }

                // Debug nodes
                if (debug_options.debug_kdtree) {
                    if (ImGui::Button("root")) {
                        debug_options.path = {};
                        if (!kdtree.empty()) {
                            debug_options.path.push_back(0);
                        }
                    }

                    if (!debug_options.path.empty() && debug_options.path.back() < kdtree.nodes().size()) {
                        auto n_index = debug_options.path.back();
                        auto const &n = kdtree.nodes().at(n_index);
                        debug_draw_node(n_index);

                        if (n.is_internal()) {
                            if (ImGui::Button("children_left")) {
                                debug_options.path.push_back(n_index + 1);
                            }
                            if (ImGui::IsItemHovered()) {
                                debug_draw_node(n_index + 1);
                            }

                            if (ImGui::Button("children_right")) {
                                debug_options.path.push_back(n.next_child());
                            }
                            if (ImGui::IsItemHovered()) {
                                debug_draw_node(n.next_child());
                            }
                        } else {
                            // Disabled version
                            ImGui::PushDisabled();
                            if (ImGui::Button("children_left")) {
                            }
                            if (ImGui::Button("children_right")) {
                            }
                            ImGui::PopDisabled();
                        }

                        if (debug_options.path.size() > 1) {
                            if (ImGui::Button("parent")) {
                                debug_options.path.pop_back();
                            }
                            if (ImGui::IsItemHovered()) {
                                debug_draw_node(debug_options.path.back());
                            }
                        } else {
                            // Disabled version
                            ImGui::PushDisabled();
                            if (ImGui::Button("parent")) {
                            }
                            ImGui::PopDisabled();
                        } {
                            // Info
                            static decltype(n_index) last_node = -1;
                            if (last_node != n_index) {
                                std::stringstream ss;
                                ss << "Node index: " << n_index << std::endl;
                                ss << "Node: " << &n << std::endl;
                                ss << "Node type: " << (n.is_internal() ? "INTERNAL" : "LEAF") << std::endl;
                                if (n.is_internal()) {
                                    ss << "Split: " << n.split() << std::endl;
                                    ss << "Axis: " << n.axis() << std::endl;
                                    ss << "\tLeft child idx: " << n_index + 1 << std::endl;
                                    ss << "\tRight child idx: " << n.next_child() << std::endl;
                                } else {
                                    ss << "Primitive start: " << n.primitive_start() << std::endl;
                                    ss << "Primitive count: " << n.primitive_count() << std::endl;
                                }
                                debug_options.node_info = ss.str();
                                last_node = n_index;
                            }
                            ImGui::PushDisabled();
                            ImGui::SetNextItemWidth(-20.0f);
                            ImGui::InputTextMultiline("a", &debug_options.node_info, ImVec2(0, 200));
                            ImGui::PopDisabled();
                        }
                    }
                }
            }
            ImGui::End();

            // Last manual ray
            if (debug_options.debug_manual_ray) {
                glEnable(GL_BLEND);
                CS350::Stats::instance().reset();
                debug.draw_segment(vp, debug_options.last_manual_ray.start,
                                   debug_options.last_manual_ray.start + debug_options.last_manual_ray.dir * 500.0f,
                                   {1, 0.25f, 1, 1});
                float t = 0.0f;
                auto const *triangle_material = scene.kdtree.get_closest(debug_options.last_manual_ray, t);
                if (triangle_material != nullptr) {
                    auto const &tri = *triangle_material;
                    debug.draw_triangle(vp, tri[0], tri[1], tri[2], {1, 0.25f, 1, 0.5f});
                }
            }
            ImGuiEndFrame();
        }
    }
    CS350::Window::destroy_system();
    return 0;
}
