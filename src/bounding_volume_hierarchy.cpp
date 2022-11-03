#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>
#include <typeinfo>

uint32_t level;
uint32_t leaves;
const uint16_t max_level = 14;

std::vector<Mesh> meshes;

uint32_t nrOfNodes;

glm::vec3 BoundingVolumeHierarchy::calculateLowerBound(std::vector<std::pair<uint32_t, uint32_t>> triangles)
{
    glm::vec3 lowerBound = glm::vec3 { std::numeric_limits<float>::max() };
    for (std::pair p : triangles) {
        glm::uvec3 triangle = meshes[p.first].triangles[p.second];
        glm::vec3 vert1 = meshes[p.first].vertices[triangle.x].position;
        glm::vec3 vert2 = meshes[p.first].vertices[triangle.y].position;
        glm::vec3 vert3 = meshes[p.first].vertices[triangle.z].position;
        lowerBound.x = std::min({ lowerBound.x, vert1.x,
                vert2.x, vert3.x });

        lowerBound.y = std::min({ lowerBound.y, vert1.y,
                vert2.y, vert3.y });

        lowerBound.z = std::min({ lowerBound.z, vert1.z,
                vert2.z, vert3.z });
    }
    
    return lowerBound;
}

glm::vec3 BoundingVolumeHierarchy::calculateUpperBound(std::vector<std::pair<uint32_t, uint32_t>> triangles)
{
    glm::vec3 upperBound = glm::vec3 { -std::numeric_limits<float>::max() };
    for (std::pair p : triangles) {
        glm::uvec3 triangle = meshes[p.first].triangles[p.second];
        glm::vec3 vert1 = meshes[p.first].vertices[triangle.x].position;
        glm::vec3 vert2 = meshes[p.first].vertices[triangle.y].position;
        glm::vec3 vert3 = meshes[p.first].vertices[triangle.z].position;
        upperBound.x = std::max({ upperBound.x, vert1.x,
            vert2.x, vert3.x });

        upperBound.y = std::max({ upperBound.y, vert1.y,
            vert2.y, vert3.y });

        upperBound.z = std::max({ upperBound.z, vert1.z,
            vert2.z, vert3.z });
    }
    return upperBound;
}

void BoundingVolumeHierarchy::createChildren(uint32_t node, uint32_t levelLocal)
{
    level = std::max( level, levelLocal );
    nodes.push_back(Node {});
    uint32_t child1 = nrOfNodes++;
    nodes.push_back(Node {});
    uint32_t child2 = nrOfNodes++;
    struct sort_by_centroid_x {
        inline bool operator()(const std::pair<uint32_t, uint32_t>& triangle1, const std::pair<uint32_t, uint32_t>& triangle2)
        {
            float centroid1 = ((meshes[triangle1.first].vertices[meshes[triangle1.first].triangles[triangle1.second].x].position + 
                meshes[triangle1.first].vertices[meshes[triangle1.first].triangles[triangle1.second].y].position + 
                meshes[triangle1.first].vertices[meshes[triangle1.first].triangles[triangle1.second].z].position) / 3.f).x;
            float centroid2 = ((meshes[triangle2.first].vertices[meshes[triangle2.first].triangles[triangle2.second].x].position + 
                meshes[triangle2.first].vertices[meshes[triangle2.first].triangles[triangle2.second].y].position + 
                meshes[triangle2.first].vertices[meshes[triangle2.first].triangles[triangle2.second].z].position) / 3.f).x;
            return (centroid1 < centroid2);
        }
    };

    struct sort_by_centroid_y {
        inline bool operator()(const std::pair<uint32_t, uint32_t>& triangle1, const std::pair<uint32_t, uint32_t>& triangle2)
        {
            float centroid1 = ((meshes[triangle1.first].vertices[meshes[triangle1.first].triangles[triangle1.second].x].position + meshes[triangle1.first].vertices[meshes[triangle1.first].triangles[triangle1.second].y].position + meshes[triangle1.first].vertices[meshes[triangle1.first].triangles[triangle1.second].z].position) / 3.f).y;
            float centroid2 = ((meshes[triangle2.first].vertices[meshes[triangle2.first].triangles[triangle2.second].x].position + meshes[triangle2.first].vertices[meshes[triangle2.first].triangles[triangle2.second].y].position + meshes[triangle2.first].vertices[meshes[triangle2.first].triangles[triangle2.second].z].position) / 3.f).y;
            return (centroid1 < centroid2);
        }
    };

    struct sort_by_centroid_z {
        inline bool operator()(const std::pair<uint32_t, uint32_t>& triangle1, const std::pair<uint32_t, uint32_t>& triangle2)
        {
            float centroid1 = ((meshes[triangle1.first].vertices[meshes[triangle1.first].triangles[triangle1.second].x].position + meshes[triangle1.first].vertices[meshes[triangle1.first].triangles[triangle1.second].y].position + meshes[triangle1.first].vertices[meshes[triangle1.first].triangles[triangle1.second].z].position) / 3.f).z;
            float centroid2 = ((meshes[triangle2.first].vertices[meshes[triangle2.first].triangles[triangle2.second].x].position + meshes[triangle2.first].vertices[meshes[triangle2.first].triangles[triangle2.second].y].position + meshes[triangle2.first].vertices[meshes[triangle2.first].triangles[triangle2.second].z].position) / 3.f).z;
            return (centroid1 < centroid2);
        }
    };
    if (levelLocal % 3 == 1) {
        std::sort(nodes[node].trianglesOrNodes.begin(), nodes[node].trianglesOrNodes.end(), sort_by_centroid_x());
    } else if (levelLocal % 3 == 2) {
        std::sort(nodes[node].trianglesOrNodes.begin(), nodes[node].trianglesOrNodes.end(), sort_by_centroid_y());
    } else {
        std::sort(nodes[node].trianglesOrNodes.begin(), nodes[node].trianglesOrNodes.end(), sort_by_centroid_z());
    }

    nodes[child1].trianglesOrNodes = std::vector<std::pair<uint32_t, uint32_t>>();
    nodes[child2].trianglesOrNodes = std::vector<std::pair<uint32_t, uint32_t>>();
    for (int i = 0; i < nodes[node].trianglesOrNodes.size() / 2; i++) {
        nodes[child1].trianglesOrNodes.push_back(nodes[node].trianglesOrNodes[i]);
    }
    for (int i = nodes[node].trianglesOrNodes.size() / 2; i < nodes[node].trianglesOrNodes.size(); i++) {
        nodes[child2].trianglesOrNodes.push_back(nodes[node].trianglesOrNodes[i]);
    }
    nodes[node].trianglesOrNodes = std::vector<std::pair<uint32_t, uint32_t>>();
    nodes[child1].lowerBound = calculateLowerBound(nodes[child1].trianglesOrNodes);
    nodes[child1].upperBound = calculateUpperBound(nodes[child1].trianglesOrNodes);

    nodes[child2].lowerBound = calculateLowerBound(nodes[child2].trianglesOrNodes);
    nodes[child2].upperBound = calculateUpperBound(nodes[child2].trianglesOrNodes);

    nodes[child1].leaf = true;
    nodes[node].trianglesOrNodes.push_back(std::pair { child1, 0 });
    nodes[child2].leaf = true;
    nodes[node].trianglesOrNodes.push_back(std::pair { child2, 0 });
    leaves += 2;

    if (levelLocal < max_level) {
        if (nodes[child1].trianglesOrNodes.size() >= 2) {
            nodes[child1].leaf = false;
            createChildren(child1, levelLocal+1);
        }

        if (nodes[child2].trianglesOrNodes.size() >= 2) {
            nodes[child2].leaf = false;
            createChildren(child2, levelLocal + 1);
        }
    }

}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
    , nodes(std::vector<Node>())
{
    level = 0;
    leaves = 0;
    nrOfNodes = 0;

    meshes = (*pScene).meshes;
    std::cout << meshes[0].triangles.size() << '\n';
    nodes.push_back(Node {});
    nrOfNodes++;
    nodes[0].lowerBound = glm::vec3 { std::numeric_limits<float>::max() };
    nodes[0].upperBound = glm::vec3 { -std::numeric_limits<float>::max() };
    for (int i = 0; i < meshes.size(); i++) {
        for (int j = 0; j < meshes[i].triangles.size(); j++) {
            nodes[0].trianglesOrNodes.push_back(std::pair { i, j });
        }
    }
    //std::cout << root.trianglesOrNodes.size() << '\n';
    nodes[0].lowerBound = calculateLowerBound(nodes[0].trianglesOrNodes);
    nodes[0].upperBound = calculateUpperBound(nodes[0].trianglesOrNodes);
    nodes.push_back(nodes[0]);
    
    if (level < max_level and nodes[0].trianglesOrNodes.size() >= 2) {
        level++;
        nodes[0].leaf = false;
        createChildren(0, 1);
    } else {
        leaves++;
    }
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return level + 1;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return leaves;
}

void BoundingVolumeHierarchy::drawLevel(int level, int index)
{
    if (level > 0 and !nodes[index].leaf) {
        drawLevel(level - 1, nodes[index].trianglesOrNodes[0].first);
        drawLevel(level - 1, nodes[index].trianglesOrNodes[1].first);
    } else if (level == 0) {
        AxisAlignedBox aabb { nodes[index].lowerBound, nodes[index].upperBound };
        drawAABB(aabb, DrawMode::Wireframe, glm::vec3(1.f, 1.0f, 1.f), 1.f);
    }
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    //AxisAlignedBox aabb { nodes[0].lowerBound, nodes[0].upperBound };
    //drawAABB(aabb, DrawMode::Wireframe);
    //drawAABB(aabb, DrawMode::Wireframe, glm::vec3(1.f, 1.0f, 1.f), 0.1f);

   drawLevel(level, 0);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    for (Node n : nodes) {
        if (n.leaf) {
            aabb = AxisAlignedBox { n.lowerBound, n.upperBound };
            if (leafIdx == 0)
                break;
            leafIdx--;
        }
    }
    
    drawAABB(aabb, DrawMode::Wireframe);
    
     //drawAABB(nodes[leafIdx].aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    //std::cout << "aabb nr: " << nodes.size() << '\n';
    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}
    // Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    bool hit = false;
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {

                    hitInfo.vertices = { v0, v1, v2 }; //store the vertices in hitInfo 
                    
                    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position,
                        v1.position, v2.position, ray.origin + ray.t * ray.direction); //update the barycentric coordinate of hitInfo
                    if (features.enableTextureMapping)
                        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord,  //update the texture coordinate of hitInfo
                            v2.texCoord, hitInfo.barycentricCoord);
                    
                    //if the normal interpolation flag enables, we update the normal of hitInfo equal to the interpolated normal. If not, we compute the cross product of vector v0->v1 and v0->v2
                    if (features.enableNormalInterp) {
                        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
                    } else {
                        hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                    }

                    hitInfo.material = mesh.material;
                    hit = true;
                }
            }
        }

        if (features.enableTextureMapping) {
            if (hitInfo.material.kdTexture) {
                hitInfo.material.kd = acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
            }
        }
        if (features.extra.enableBilinearTextureFiltering) {
            if (hitInfo.material.kdTexture) {
                hitInfo.material.kd = acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
            }
        }
        
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit || recursiveTraversal(0, ray, hitInfo, features);
    }   

    
}

bool BoundingVolumeHierarchy::recursiveTraversal(int nodeIndex, Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    bool hit = false;
    if (nodes[nodeIndex].leaf) {
        for (std::pair<uint32_t, uint32_t> triangle : nodes[nodeIndex].trianglesOrNodes) {
            glm::uvec3 tri = meshes[triangle.first].triangles[triangle.second];
            const auto v0 = meshes[triangle.first].vertices[tri[0]];
            const auto v1 = meshes[triangle.first].vertices[tri[1]];
            const auto v2 = meshes[triangle.first].vertices[tri[2]];

            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {

                hitInfo.vertices = { v0, v1, v2 }; // store the vertices in hitInfo

                hitInfo.barycentricCoord = computeBarycentricCoord(v0.position,
                    v1.position, v2.position, ray.origin + ray.t * ray.direction); // update the barycentric coordinate of hitInfo
                if (features.enableTextureMapping)
                    hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, // update the texture coordinate of hitInfo
                        v2.texCoord, hitInfo.barycentricCoord);

                // if the normal interpolation flag enables, we update the normal of hitInfo equal to the interpolated normal. If not, we compute the cross product of vector v0->v1 and v0->v2
                if (features.enableNormalInterp) {
                    hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
                } else {
                    hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                }

                hitInfo.material = meshes[triangle.first].material;
                hit = true;
            }
        }
        return hit;
    } else {
        AxisAlignedBox aabbNode = AxisAlignedBox { nodes[nodeIndex].lowerBound, nodes[nodeIndex].upperBound };
        float tBackup = ray.t;
        if (intersectRayWithShape(aabbNode, ray) ){
            drawRay(ray);
            ray.t = tBackup;
            hit |= recursiveTraversal(nodes[nodeIndex].trianglesOrNodes[0].first, ray, hitInfo, features);
            hit|= recursiveTraversal(nodes[nodeIndex].trianglesOrNodes[1].first, ray, hitInfo, features);
            if (!hit)
                ray.t = tBackup;
            return hit;
        } else {
            ray.t = tBackup;
            return false;
        }
    }
}