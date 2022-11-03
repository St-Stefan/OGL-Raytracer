#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#include "texture.h"
#ifdef NDEBUG
#include <omp.h>
#endif

Image negx = Image("../../../data/negx.jpg");
Image negy = Image("../../../data/negy.jpg");
Image negz = Image("../../../data/negz.jpg");
Image posx = Image("../../../data/posx.jpg");
Image posy = Image("../../../data/posy.jpg");
Image posz = Image("../../../data/posz.jpg");

std::vector<Image> image = { posx, negx, posy, negy, posz, negz };

//Implementing the recursive ray-tracer 
//Source: Chapter 4.8, Fundamentals of Computer Graphics, Fourth Edition.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        Ray reflection = computeReflectionRay(ray, hitInfo);
        if (features.extra.enableEnvironmentMapping) {
            hitInfo.material.kd = environmentMapping(image, reflection, features);
        }

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            if (rayDepth < 10 && hitInfo.material.ks != glm::vec3{ 0.0f }) {
                Lo += hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
            }
        }
        //if normal interpolation flag enables, we draw normal at 3 vertices and the interpolated normal.    
        /* if (features.enableNormalInterp) {
            drawRay(Ray { hitInfo.vertices[0].position, hitInfo.vertices[0].normal, 2 }, glm::vec3 { 1.0f, 0.0f, 0.0f });
            drawRay(Ray { hitInfo.vertices[1].position, hitInfo.vertices[1].normal, 2 }, glm::vec3 { 1.0f, 0.0f, 0.0f });
            drawRay(Ray { hitInfo.vertices[2].position, hitInfo.vertices[2].normal, 2 }, glm::vec3 { 1.0f, 0.0f, 0.0f });
            drawRay(Ray { ray.origin + ray.t * ray.direction,
                      hitInfo.normal, 2 },
             glm::vec3 { 0.0f, 1.0f, 0.0f });
        }*/

        // Visual Debug: Draw a ray with a color which is the returned value from computeLightContribution
        drawRay(ray, Lo);
        // Set the color of the pixel to white if the ray hits.
        return Lo;

    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        if (features.extra.enableEnvironmentMapping) {

            return environmentMapping(image, ray, features);
        } else {
            return glm::vec3(0.0f);
        }
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
        }
    }
}