#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

//Compute the Phong Shading Model 
//Source: https://users.cs.northwestern.edu/~ago820/cs395/Papers/Phong_1975.pdf
const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{

    if (features.enableShading) {
        // Compute direction of reflected vector
        glm::vec3 intersectionPoint = ray.origin + ray.direction * ray.t;
        glm::vec3 directionOfIncomingRay = glm::normalize(lightPosition - intersectionPoint);
        glm::vec3 normal = glm::normalize(hitInfo.normal);
        glm::vec3 rayFromEyeDirectionNormalized = -glm::normalize(ray.direction);
        glm::vec3 reflectedVector = glm::normalize(2 * (glm::dot(normal, directionOfIncomingRay)) * normal - directionOfIncomingRay);
        
        glm::vec3 diffuseTerm;
        // Diffuse Term
        if (glm::dot(normal, directionOfIncomingRay) < 0) {
            diffuseTerm = glm::vec3 { 0.0f };
        } else {
            if (features.enableTextureMapping) {
                diffuseTerm = lightColor
                    * acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features) // Based on the book Fundamentals of Computer Graphics, chapter 11.1
                                                                                                  //, we will replace the value Kd = acquireTexel(...)  
                    * glm::dot(normal, directionOfIncomingRay);
            } else {
                diffuseTerm = lightColor
                    * hitInfo.material.kd
                    * glm::dot(normal, directionOfIncomingRay);
            }
        }

        glm::vec3 specularTerm;
        // Specular Term
        if (glm::dot(rayFromEyeDirectionNormalized, reflectedVector) < 0) {
            specularTerm = glm::vec3 { 0.0f };
        } else {
            specularTerm = lightColor
                * hitInfo.material.ks
                * glm::pow(glm::dot(rayFromEyeDirectionNormalized, reflectedVector), hitInfo.material.shininess);
        }

        return diffuseTerm + specularTerm;
    } else {
        if (features.enableShading) {
            return acquireTexel(*hitInfo.material.kdTexture.get(), hitInfo.texCoord, features);
        }
        return hitInfo.material.kd;
    }
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    glm::vec3 normal = glm::normalize(hitInfo.normal);
    glm::vec3 direction = -glm::normalize(ray.direction);
    glm::vec3 intersectionPoint = ray.origin + ray.direction * ray.t;
    glm::vec3 reflectionDirection = glm::normalize(-direction + 2 * (glm::dot(direction, normal)) * normal);

    
    // TODO: implement the reflection ray computation.
    //Add the offset to the intersectionPoint
    return Ray {intersectionPoint + reflectionDirection * 0.001f, reflectionDirection};
}