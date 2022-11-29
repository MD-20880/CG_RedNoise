#include "RayTriangleIntersection.h"

RayTriangleIntersection::RayTriangleIntersection() = default;
RayTriangleIntersection::RayTriangleIntersection(const glm::vec3 &point, float distance, const ModelTriangle &triangle, size_t index) :
		intersectionPoint(point),
		distanceFromCamera(distance),
		intersectedTriangle(triangle),
		triangleIndex(index),
		possibleSolution(glm::vec3(-1,-1,-1)){}
RayTriangleIntersection::RayTriangleIntersection(const glm::vec3 &point, float distance, const ModelTriangle &triangle, size_t index,const glm::vec3 &possibleSolution):
		intersectionPoint(point),
		distanceFromCamera(distance),
		intersectedTriangle(triangle),
		triangleIndex(index),
		possibleSolution(possibleSolution){}


std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection) {
	os << "Intersection is at [" << intersection.intersectionPoint[0] << "," << intersection.intersectionPoint[1] << "," <<
	   intersection.intersectionPoint[2] << "] on triangle " << intersection.intersectedTriangle <<
	   " at a distance of " << intersection.distanceFromCamera;
	return os;
}
