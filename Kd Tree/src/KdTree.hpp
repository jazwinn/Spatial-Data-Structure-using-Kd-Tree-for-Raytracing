/**
 * @file
 *  KdTree.hpp
 * @author
 *  Jaz Winn Ng, 670001224, jazwinn.ng@digipen.edu
 * @date
 *  2025/05/31
 * @brief
 *  Provides the declaration and definition of the KdTree and
 *  Ray query
 * @copyright
 *  Copyright (C) 2025 DigiPen Institute of Technology.
 */

#ifndef KDTREE_HPP
#define KDTREE_HPP

#include "Math.hpp"
#include "Shapes.hpp"
#include "Stats.hpp"
#include "Logging.hpp"

#include <glm/gtc/epsilon.hpp>
#include <vector>
#include <iostream>
#include <functional>
#include <sstream>
#include <stack>

namespace CS350 {

    constexpr float epsilon = 1e-5f;

    /**
     * Construction configuration
     */
    struct KdTreeConfig {
        float cost_traversal    = 1.0f;
        float cost_intersection = 80.0f;
        int   max_depth         = 5;
        int   min_triangles     = 50;
        int   max_split         = 100;
    };

    /**
     * KDTree Node structure
     *  Internals: Split, Axis, Right subnode index. (left subnode is next in the array)
     *  Leafs: Index of first triangle, amount of triangles
     */
    struct KdTreeNode {

    public:
        /**
         * @brief
		 *  Returns whether the node is an internal node or a leaf node.
         * @return
		 *  Returns true if the node is an internal node, false if it is a leaf node.
         */
        bool is_internal() const {
            return (m_count & 0b11) != 0b11;
        }

        /**
         * @brief
         *  Returns whether the node is an internal node or a leaf node.
         * @return
         *  Returns true if the node is an leaf node, false if it is a internal node.
         */
        bool is_leaf() const {
            return (m_count & 0b11) == 0b11;
        }

		/**
		 * @brief
		 *  Returns the index of the right sub node.
		 * @return
		 *  Returns the index of the right sub node.
		 */
        size_t next_child() const {

            return m_subnode_index >> 2;
        }

		/**
		 * @brief
		 *  Returns the number of primitives in this node.
		 * @return
		 *  Number of primitives in this node.
		 */
        int primitive_count() const {
            return m_count >> 2;
        }

		/**
		 * @brief
		 *  Returns the index of the first primitive in this node.
		 * @return
		 *  Index of the first primitive in this node.
		 */
        unsigned primitive_start() const {
            return m_start_primitive;
        }

        /**
         * @brief
         *  Returns the axis of the split point.
		 *  x = 0, y = 1, z = 2
         * @return
         *  axis of the split point.
         */
        size_t axis() const {
            return static_cast<size_t>(m_subnode_index & 0b11);
        }

        /**
         * @brief
         *  Returns the split point of an axis.
         * @return
         *  Split point of an axis.
         */
		float split() const {
			return m_split;
		}

		/**
		 * @brief
		 *  Sets the split point of an axis.
		 * @param split
		 *  Split point of an axis.
		 */
        void Set_Split(float split) {
            m_split = split;
        }

		/**
		 * @brief
		 *  Sets the index of the first primitive in this node.
		 * @param start
		 *  Index of the first primitive in this node.
		 */
        void Set_Start_Primitive(unsigned start) {
            m_start_primitive = start;
        }

        /**
         * @brief
         *  Sets the index of the right child.
         * @param index
         *  Index of the right child.
         */
        void Set_Subnode_Index(unsigned index) {
            unsigned axis = m_subnode_index & 0b11;
            m_subnode_index = (index << 2) | axis;
        }

        /**
         * @brief
		 *  Sets the axis of the split point.
         * @param axis
		 *  Axis of the split point.
         */
		void Set_Axis(unsigned axis) {

			unsigned index = m_subnode_index & (~0b11);
			m_subnode_index = index | (axis & 0b11);
		}

		/**
		 * @brief
		 *  Sets whether the node is a leaf or not.
		 * @param leaf
		 *  True if the node is a leaf, false if it is an internal node.
		 */
		void Set_Leaf(bool leaf) {
            if (leaf) {
				m_subnode_index = m_subnode_index | 0b11;
            }
            else {
                m_subnode_index = m_subnode_index & (~0b11);
            }

		}

		/**
		 * @brief
		 *  Set the number of triangles in the node.
		 * @param count
		 *  Number of triangles in the node.
		 */
		void Set_Count(unsigned count) {
			unsigned leaf = m_count & 0b11;
            m_count = (count << 2) | leaf;
		}

      private:   

		  union
		  {
			  float m_split;              // Split position (internal only)
			  unsigned m_start_primitive; // Index of first triangle (leafs only)
		  };

		  union
		  {
			  unsigned m_subnode_index;   // Subnode index (30MSB) + Axis (2LSB) (internal only)
			  unsigned m_count;           // Amount of triangles (30MSB) + 0b11 (leafs only)
		  };
          
    };

    template <typename T>
    class KdTree {
        std::vector<T>          m_triangles; // All recorded triangles (may contain duplicates)
        std::vector<KdTreeNode> m_nodes;     // KDTree nodes
        std::vector<Aabb>       m_aabbs;     // AABBs of nodes (same order)
        KdTreeConfig            m_cfg;       // Configuration
      public:

		/**
		* @brief
		*  Constructs a KdTree with the given configuration.
		* @param all_triangles
		*  The triangles to build the KdTree with.
		*/
        void build(std::vector<T> const& all_triangles);

		/**
		 * @brief
		 *  Query Ray against the KdTree to find closest triangle intersection.
		 * @param r
		 *  Ray of the triangle.
		 * @param out_t
		 *  Closest intersection distance from the ray start.
		 * @return
		 *  First Triangle that intersects the ray.
		 */
        T const* get_closest(Ray const& r, float& out_t) const;

        std::ostream& dump(std::ostream&) const;
        std::ostream& dump_graph(std::ostream&) const;

		// Accessors
        auto& nodes() const { return m_nodes; }
        auto& aabbs() const { return m_aabbs; }
        bool  empty() const { return m_triangles.empty(); }
        auto& cfg() const { return m_cfg; }
        auto& cfg() { return m_cfg; }
        auto& triangles() const { return m_triangles; }
        auto& triangles() { return m_triangles; } // Debug

		/**
		 * @brief
		 *  Returns the depth of the KdTree.
		 * @return
		 *  Depth of the KdTree, -1 if empty.
		 */
        int   get_depth() const { if (m_nodes.empty() ) return -1; return get_depth_rec(0); }

        /**
         * @brief
		 *  Gets all triangles from node index down to its leafs.
		 * @param out_triangles
		 *  All the triangles in the leafs below the node index.
         */
        void  get_triangles(size_t node_index, std::vector<T>& out_triangles) const;


      private:

        /**
        * @brief
        *  Builds the KdTree recursively, called by the build function.
        * @param all_triangles
        *  The triangles to build the KdTree with.
        * @param parentIndex
        *  The index of the parent node.
        * @param level
        *  The current level of the tree.
        */
        void Build_Recur(std::vector<T> const& all_triangles, size_t currentNodeIndex, int level = 1);

        /**
         * @brief
         *  Finds the depth of the KdTree recursively, called by the get_depth function.
         * @param node_index
         *  The index of the node.
         */
        int get_depth_rec(int node_index) const;
      
    };

    template <typename T>
    std::ostream& KdTree<T>::dump(std::ostream& os) const {
        std::function<void(int n, int level)> node_dump;
        node_dump = [&](int n, int level) {
            auto const& node = m_nodes.at(n);
            std::string tab(level * 2, ' ');
            os << tab << "Node " << n << "[";
            if (node.is_internal()) {
                os << "internal, "
                   << "SplitPosition at " << char('x' + node.axis())
                   << "=" << node.split()
                   << "]\n";
                node_dump(n + 1, level + 1);
                node_dump(static_cast<int>(node.next_child()), level + 1);
            } else {
                os << "leaf, "
                   << node.primitive_start() << ":" << node.primitive_start() + node.primitive_count()
                   << "]\n";
            }
        };
        node_dump(0, 0);
        return os;
    }

    template <typename T>
    std::ostream& KdTree<T>::dump_graph(std::ostream& os) const {
        std::function<int(int n)> node_triangle_count;
        node_triangle_count = [&](int n) {
            auto const& node = m_nodes.at(n);
            if (node.is_internal()) {
                return node_triangle_count(n + 1) + node_triangle_count(static_cast<int>(node.next_child()));
            } else {
                return node.primitive_count();
            }
        };

        std::function<void(int n, int n_parent)> node_dump;
        node_dump = [&](int n, int n_parent) {
            std::stringstream s;
            s << "NODE" << n;
            std::string name = s.str();
            os << "\t" << name << "[label=\"";
            auto const& node = m_nodes.at(n);
            if (node.is_internal()) {
                os << "SplitPosition: " << char('x' + node.axis()) << " at " << node.split() << "\\n"
                   << node_triangle_count(n) << " subtriangles";
            } else {
                os << node.primitive_count() << " triangles";
            }
            os << "\"];\n";

            // Edges
            if (n) {
                s = std::stringstream();
                s << "NODE" << n_parent;
                std::string parent_name = s.str();
                os << "\t" << parent_name << " -> " << name << ";\n";
            }

            // Children
            if (node.is_internal()) {
                node_dump(n + 1, n);
                node_dump(static_cast<int>(node.next_child()), n);
            }
        };
        os << "digraph kdtree {\n";
        os << "\tnode[group=\"\", shape=none, style=\"rounded,filled\", fontcolor=\"#101010\"]\n";
        node_dump(0, 0);
        os << "}";
        return os;
    }

    template<typename T>
	int KdTree<T>::get_depth_rec(int node_index) const {
		
		//recurse down the tree until we reach a leaf node
		if (m_nodes[node_index].is_leaf()) {
			return 0;
		}


        int indexChld0 = node_index + 1;
        int indexChld1 = static_cast<int>(m_nodes[node_index].next_child());


        int child0Depth = get_depth_rec(indexChld0);

		int child1Depth = get_depth_rec(indexChld1);


		return 1 + glm::max(child0Depth, child1Depth);


	}



	template<typename T>
	void KdTree<T>::build(std::vector<T> const& all_triangles) {

        if (all_triangles.empty()) {
            return;
        }
        m_triangles.clear();
        m_nodes.clear();
        m_aabbs.clear();
        
        //create initial node
        m_nodes.emplace_back(KdTreeNode{});

		//create AABB
		vec3 min = all_triangles[0][0];
		vec3 max = all_triangles[0][0];
		for (size_t n{}; n < all_triangles.size(); n++) {

			vec3 minPoint = glm::min(all_triangles[n][0], all_triangles[n][1], all_triangles[n][2]);
			vec3 maxPoint = glm::max(all_triangles[n][0], all_triangles[n][1], all_triangles[n][2]);


			min = glm::min(min, minPoint);
			max = glm::max(max, maxPoint);
		}
		m_aabbs.emplace_back(min, max);
		
        Build_Recur(all_triangles, 0);

	}


    template<typename T>
    void KdTree<T>::Build_Recur(std::vector<T> const& all_triangles, size_t parentIndex, int level) {



        //test config
        // if condition is true, turn node to leaf
        if (level >= m_cfg.max_depth || all_triangles.size() <= static_cast<size_t>(m_cfg.min_triangles)) {
            size_t index = m_triangles.size();
            m_triangles.insert(m_triangles.end(), all_triangles.begin(), all_triangles.end());


            m_nodes[parentIndex].Set_Start_Primitive(static_cast<unsigned int>(index));
            m_nodes[parentIndex].Set_Count(static_cast<unsigned int>(all_triangles.size()));
            m_nodes[parentIndex].Set_Leaf(true);

            return;
        }

        Aabb* currentAabb = &m_aabbs.back();

        //get best axis
        int bestAxis = currentAabb->longest_axis();


        //conduct SAH
        //Best scenario
        float lowestCost = std::numeric_limits<float>::max();
        float bestPoint = 0;

         //Store the split points (Min/Max bound of triangle)
        //Pair first being the min, second the max
        std::vector<std::pair<float, float>> triBound;
        triBound.reserve(all_triangles.size());

		//Store split points
        //Pair first being split point , second being a boolean, 0 = min, 1 = max, of the triangle
		std::vector<std::pair<float, bool>> splitPoints;
		splitPoints.reserve(all_triangles.size() * 2);


        // get all split points of the axis
		for (const auto& triangle : all_triangles) {
            triBound.push_back({ glm::min(triangle[0][bestAxis], triangle[1][bestAxis], triangle[2][bestAxis]), glm::max(triangle[0][bestAxis], triangle[1][bestAxis], triangle[2][bestAxis])});

			splitPoints.push_back({ triBound.back().first, 0});
			splitPoints.push_back({ triBound.back().second, 1});

		}

		auto sortSplit = [](const auto& lhs, const auto& rhs) {

			if (lhs.first == rhs.first) {
				return lhs.second < rhs.second;
			}
			return lhs.first < rhs.first;
			};

		std::sort(splitPoints.begin(), splitPoints.end(), sortSplit);


		//static calculation
		float parentNodeSA = currentAabb->surface_area();


        //Start with 0 object on the left, all of it on the right
		size_t leftSide = 0;
		size_t rightSide = all_triangles.size();


		// test all split points
		for (size_t n = 0; n < splitPoints.size(); n++) {

            float point = splitPoints[n].first;
			bool startOrEnd = splitPoints[n].second;
			
            // increment if split line "touches" the min triangle boundary
			if (startOrEnd == 0) {
				leftSide++;
			}


			//skip points that may be outside/on current aabb
			if (point <= currentAabb->min[bestAxis] || point >= currentAabb->max[bestAxis]) {
				continue;
			}

            //Perform cost heuristic
			vec3 leftMax = currentAabb->max;
			leftMax[bestAxis] = point;

			Aabb leftAabb{ currentAabb->min, leftMax };

			vec3 rightMin = currentAabb->min;
			rightMin[bestAxis] = point;

			Aabb rightAabb{ rightMin, currentAabb->max };

			float probLeft = leftAabb.surface_area() / parentNodeSA;
			float probRight = rightAabb.surface_area() / parentNodeSA;

			float cost = m_cfg.cost_traversal + m_cfg.cost_intersection * (probLeft * static_cast<float>(leftSide) + probRight * static_cast<float>(rightSide));

            // decrement max triangle boundary "leaves" the split line
			if (startOrEnd == 1) {
				rightSide--;
			}


			//find intersecting triangles
		    size_t intersect = (leftSide + rightSide) - all_triangles.size();
            size_t leftOnly = leftSide - intersect;
            size_t rightOnly = rightSide - intersect;

            // Check that each side have at least 20 unique triang;es
			if (leftOnly < static_cast<size_t>(m_cfg.min_triangles) || rightOnly < static_cast<size_t>(m_cfg.min_triangles)) {
				continue;
			}

            //keep track of lowest cost
			if (cost < lowestCost) {
				lowestCost = cost;
				bestPoint = point;
			}


		}



        //test possibility of creating a leaf node
        float newLeafCost = m_cfg.cost_intersection * static_cast<float>(all_triangles.size());

        if (newLeafCost < lowestCost) {
            //create leaf node
            size_t index = m_triangles.size();
            m_triangles.insert(m_triangles.end(), all_triangles.begin(), all_triangles.end());


            m_nodes[parentIndex].Set_Start_Primitive(static_cast<unsigned int>(index));
            m_nodes[parentIndex].Set_Count(static_cast<unsigned int>(all_triangles.size()));
            m_nodes[parentIndex].Set_Leaf(true);


        }
        else {

			

            //generate left and right vector
            std::vector<T> leftTriangles;
            leftTriangles.reserve(all_triangles.size()/2);
            std::vector<T> rightTriangles;
            rightTriangles.reserve(all_triangles.size()/2);

			for (size_t i = 0; i < all_triangles.size(); ++i) {
				const auto& tri = all_triangles[i];

				float triMin = triBound[i].first;
				float triMax = triBound[i].second;

				if (triMin < bestPoint) {
                    leftTriangles.push_back(tri);
				}

				if (triMax > bestPoint) {
                    rightTriangles.push_back(tri);
				}

			}


            //Calculate the new AABB
			vec3 leftMax = currentAabb->max;
			leftMax[bestAxis] = bestPoint;
			Aabb bestLeftAabb = { currentAabb->min, leftMax };

			vec3 rightMin = currentAabb->min;
			rightMin[bestAxis] = bestPoint;
			Aabb bestRightAabb = { rightMin, currentAabb->max };


            size_t child0Index = m_nodes.size();
            m_nodes.emplace_back(KdTreeNode{});
            m_aabbs.emplace_back(bestLeftAabb);
            Build_Recur(leftTriangles, child0Index, level + 1);

            size_t child1Index = m_nodes.size();
            m_nodes.emplace_back(KdTreeNode{});
            m_aabbs.emplace_back(bestRightAabb);
            Build_Recur(rightTriangles, child1Index, level + 1);

            m_nodes[parentIndex].Set_Split(bestPoint);
            m_nodes[parentIndex].Set_Subnode_Index(static_cast<unsigned int>(child1Index));
            m_nodes[parentIndex].Set_Axis(bestAxis);


        }
    }

	template<typename T>
	T const* KdTree<T>::get_closest(Ray const& r, float& out_t) const {



        std::stack<size_t> stack;
        stack.push(0);

		size_t bestTriangleIndex = 0;
		float smallestIntersect = std::numeric_limits<float>::max();
        while (!stack.empty()) {
            Stats::instance().query_count++;

            size_t nodeIndex = stack.top();
            stack.pop();

            //Ray intersect and reject nodes
			Ray::Intersection intersectResult = r.intersect(m_aabbs[nodeIndex]);
			if (!intersectResult || intersectResult.tMax < epsilon ||
                intersectResult.t > smallestIntersect) {
				continue;
			}

            //Begin Traversal 
            Stats::instance().kdtree_last_traversed_nodes.push_back(nodeIndex);

			const KdTreeNode* currentNode = &m_nodes[nodeIndex];
			size_t axis = currentNode->axis();

			float tMin = intersectResult.t;
			float tMax = intersectResult.tMax;

            //If leaf is reached, check triangles
            if (currentNode->is_leaf()) {

                int count = currentNode->primitive_count();
                bool hasIntersected = false;
                for (int n{}; n < count; n++) {
                    
                    size_t index = currentNode->primitive_start() + n;
                    const T* triangle = &m_triangles[index];

                    Ray::Intersection intersect = r.intersect(static_cast<const Triangle>(*triangle));

					//end checking other nodes once we found a triangle that intersects
                    if (intersect) {

                        if (intersect.t < smallestIntersect) {
                            bestTriangleIndex = index;
                            smallestIntersect = intersect.t;
                            hasIntersected = true;
                        }
                    }

                }
                if (hasIntersected) {
                    break;
                }

                continue;
            }


             

			//check if dir is parallel to split line
			size_t leftChild = nodeIndex + 1;
			size_t rightChild = currentNode->next_child();

            float rayDir = r.dir[static_cast<int>(axis)];

			bool isRayNeg = rayDir < 0.0f;

			size_t nearIndex = isRayNeg ? rightChild : leftChild;
			size_t farIndex = isRayNeg ? leftChild : rightChild;

            //if ray is parallel, check near 
			if (glm::abs(rayDir) < epsilon) {

				stack.push(nearIndex);  
                continue;
			}

            float tSplit = (currentNode->split() - r.start[static_cast<int>(axis)]) / rayDir;
            
			if (tSplit < tMin) {
                // Ray passes split point outside of Aabb, before intersecting
				stack.push(farIndex);
			}
			else if (tSplit > tMax) {
                // Ray pass split point outside of Aabb, after intersecting
				stack.push(nearIndex);
			}
			else {
				// Ray intersects split plane
				stack.push(farIndex);
				stack.push(nearIndex);
			}
            
        }
        

        //check if ray have any intersection before returning
		if (smallestIntersect < std::numeric_limits<float>::max()) {
			out_t = smallestIntersect;
			return &m_triangles[bestTriangleIndex];
		}
		else {
            return nullptr;
		}
		
	}

	template<typename T>
	void  KdTree<T>::get_triangles(size_t node_index, std::vector<T>& out_triangles) const {
		const KdTreeNode* node = &m_nodes[node_index];

        if (node->is_leaf()) {
			out_triangles.insert(out_triangles.end(), m_triangles.begin() + node->primitive_start(), m_triangles.begin() + node->primitive_start() + node->primitive_count());
        }
        else {
            get_triangles(node_index + 1, out_triangles);
            get_triangles(node->next_child(), out_triangles);
        }
	}
}

#endif // KDTREE_HPP
