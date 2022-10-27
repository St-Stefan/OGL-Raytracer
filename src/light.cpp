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
    glm::vec3 point1 = segmentLight.endpoint0;
    glm::vec3 point2 = segmentLight.endpoint1;

    glm::vec3 lightVector = point2 - point1;

    float length = glm::length(lightVector);
    float sampleJitter = (float)rand() / (float)(RAND_MAX / length);
    
    glm::vec3 samplepoint = point1 + glm::normalize(lightVector) * sampleJitter;

    Ray tes;
    tes.origin = point1;
    tes.direction = glm::normalize(lightVector);
    tes.t = sampleJitter;
    drawRay(tes,glm::vec3(1,1,0));
    
    

    position = samplepoint;
    color = glm::vec3(0.0);
    if (point1 != point2)
        color = (length - sampleJitter) / length * segmentLight.color0 + (sampleJitter) / length * segmentLight.color1;
    // TODO: implement this function.
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    float length1 = glm::length(parallelogramLight.edge01);
    float sampleJitter1 = (float)rand() / (float)(RAND_MAX / length1);

    float length2 = glm::length(parallelogramLight.edge02);
    float sampleJitter2 = (float)rand() / (float)(RAND_MAX / length2);

    glm::vec3 samplePoint = parallelogramLight.v0 + glm::normalize(parallelogramLight.edge01) * sampleJitter1 + glm::normalize(parallelogramLight.edge02)*sampleJitter2;

    position = samplePoint;
    float coef1 = ((length1 - sampleJitter1) * (length2 - sampleJitter2)) / (length1 * length2);
    float coef2 = ((length1 - sampleJitter1) * (sampleJitter2)) / (length1 * length2);
    float coef3 = ((sampleJitter1) * (length2 - sampleJitter2)) / (length1 * length2);
    float coef4 = ((sampleJitter1) * (sampleJitter2)) / (length1 * length2);
   

    color = parallelogramLight.color0*coef1+parallelogramLight.color1*coef3+parallelogramLight.color2*coef2+parallelogramLight.color3*coef4;
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
        drawRay(lightTest, glm::vec3(0, 0, 1)); // show shadow ray in blue
        return 0.0;
    }
    else {
        drawRay(lightTest, debugColor); // otherwise shade ray the same color as the light
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
    int samplecount=250;
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
            } else if (std::holds_alternative<SegmentLight>(light)) {
                glm::vec3 segmentLightContribution { 0.0f };
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                for (int i = 0; i < samplecount; i++) {
                    glm::vec3 positi = { 0, 0, 0 };
                    glm::vec3 color = { 0, 0, 0 };
                    sampleSegmentLight(segmentLight, positi, color);
                    if (testVisibilityLightSample(positi, glm::vec3(0, 0, 1), bvh, features, ray, hitInfo) && features.enableSoftShadow || !features.enableSoftShadow) {
                        segmentLightContribution += computeShading(positi, color, features, ray, hitInfo);
                    }
                    }
                segmentLightContribution = segmentLightContribution / float(samplecount);
                lightContribution += segmentLightContribution;

                } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                glm::vec3 pLightContribution { 0, 0, 0 };
                for (int i = 0; i < samplecount; i++) {
                    glm::vec3 positi = { 0, 0, 0 };
                    glm::vec3 color = { 0, 0, 0 };
                    sampleParallelogramLight(parallelogramLight, positi, color);
                    if (testVisibilityLightSample(positi, color, bvh, features, ray, hitInfo) && features.enableSoftShadow || !features.enableSoftShadow) {
                        pLightContribution += computeShading(positi, color, features, ray, hitInfo);
                    }
                }
                pLightContribution = pLightContribution / float(samplecount);
                lightContribution += pLightContribution;
            }
        }
        return lightContribution;

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
