#include "interpolate.h"
#include <glm/geometric.hpp>

//Compute barycentric coordinates 
//Source: Chapter 2.7.2, Fundamentals of Computer Graphics, Fourth Edition.
glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // TODO: implement this function.
    glm::vec3 n = glm::cross(v1 - v0, v2 - v0); // n = (b-a) x (c-a)
    glm::vec3 n_a = glm::cross(v2 - v1, p - v1);
    glm::vec3 n_b = glm::cross(v0 - v2, p - v2);
    glm::vec3 n_c = glm::cross(v1 - v0, p - v0);

    float alpha = glm::dot(n, n_a) / glm::pow(glm::length(n), 2);
    float beta = glm::dot(n, n_b) / glm::pow(glm::length(n), 2);
    float gamma = glm::dot(n, n_c) / glm::pow(glm::length(n), 2);

    return glm::vec3(alpha, beta, gamma);
}

//Compute the interpolated normal 
glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    // TODO: implement this function.
    return barycentricCoord.x * n0 + barycentricCoord.y * n1 + barycentricCoord.z * n2;
}

//Compute the interpolated texture coordinate
glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    // TODO: implement this function.
    return barycentricCoord.x * t0 + barycentricCoord.y * t1 + barycentricCoord.z * t2;
}
