#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;



struct Node {
    glm::vec3 lowerBound, upperBound;
    bool leaf = true;
    std::vector<std::pair<uint32_t, uint32_t>> trianglesOrNodes;
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    bool recursiveTraversal(int nodeIndex, Ray& ray, HitInfo& hitInfo, const Features& features) const;

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;
    void createChildren(uint32_t node, uint32_t levelLocal);
    glm::vec3 calculateUpperBound(std::vector<std::pair<uint32_t, uint32_t>> triangles);
    glm::vec3 calculateLowerBound(std::vector<std::pair<uint32_t, uint32_t>> triangles);
    void drawLevel(int level, int index);

private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<Node> nodes;
};