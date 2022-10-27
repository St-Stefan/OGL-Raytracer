#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>


// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    position = glm::vec3(0.0);
    color = glm::vec3(0.0);
    // TODO: implement this function.
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    position = glm::vec3(0.0);
    color = glm::vec3(0.0);
    // TODO: implement this function.
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{

    float offset = 0.001; //removes shadow acne
    glm::vec3 intersectionPoint = ray.origin + ray.direction * ray.t;
    glm::vec3 rayDir = samplePos - intersectionPoint;
    Ray lightTest {glm::normalize(rayDir),intersectionPoint};
    lightTest.direction = glm::normalize(rayDir);
    lightTest.origin = intersectionPoint+offset*lightTest.direction;
    lightTest.t = glm::length(rayDir);
    HitInfo lightHitInfo;
    float shadow = false;

    if (dot(hitInfo.normal, lightTest.direction) * dot(hitInfo.normal, -ray.direction)<0)   //check if view ray and shadow ray are on the same side
        shadow = true;

    if (bvh.intersect(lightTest, lightHitInfo, features)) {     //if there is an intersection regard point as in shadow
        if (lightTest.t < glm::length(rayDir) - offset)         //subtract the offset from the final length
        shadow = true;
    }

    if (shadow) {
        drawRay(lightTest, debugColor);     //if in shadow ray is shown in chosen colour
        return 0.0;
    }
    else {
        drawRay(lightTest, glm::vec3(0, 0, 1)); // otherwise shade as yellow
            return 1.0;
    }
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 lightContribution {0.0f};
    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.
        for (const auto& light : scene.lights) {
            //only calculate shadows when flag is set to true
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                if (testVisibilityLightSample(pointLight.position, glm::vec3(0, 0, 1), bvh, features, ray, hitInfo) && features.enableHardShadow || !features.enableHardShadow) {
                    lightContribution += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                    }
            }
        }
        return lightContribution;

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
