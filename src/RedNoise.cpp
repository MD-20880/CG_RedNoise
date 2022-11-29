#include <CanvasTriangle.h>
#include "ModelTriangle.h"
#include <CanvasPoint.h>
#include <DrawingWindow.h>
#include <TextureMap.h>
#include <Colour.h>
#include <Utils.h>
#include <thread>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>	
#include <iostream>
#include <fstream>
#include <map>
#include "RayTriangleIntersection.h"
#include <glm/gtx/string_cast.hpp>

#define WIDTH 320
#define HEIGHT 240
#define PI 3.14

uint32_t fromColour(Colour colour){

	uint32_t colour_in_uint32 = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
	return colour_in_uint32;
}


std::vector<std::vector<float>> makedBuffer(){
	std::vector<std::vector<float>> result;
	for (size_t i = 0; i < HEIGHT; i++)
	{
		std::vector<float> row;
		for (size_t j = 0; j < WIDTH; j++)
		{
			row.push_back(FLT_MIN);
		}
		result.push_back(row);
	}
	return result;
}

std::vector<std::vector<uint32_t>> makepBuffer(){
	std::vector<std::vector<uint32_t>> result;
	for (size_t i = 0; i < HEIGHT; i++)
	{
		std::vector<uint32_t> row;
		for (size_t j = 0; j < WIDTH; j++)
		{
			row.push_back(fromColour(Colour(0,0,0)));
		}
		result.push_back(row);
	}
	return result;
}

//Depth Buffer of each pixel
std::vector<std::vector<float>> dbuffer = makedBuffer();
//PixelBuffer that Could be Updated before running
std::vector<std::vector<uint32_t>> pbuffer = makepBuffer();
glm::vec3 campos(0,0,4);
uint8_t state = 0x01;


void calculateNormal(std::vector<ModelTriangle> &list){
	for (int i = 0; i <list.size(); i++){
		ModelTriangle triangle = list[i];
		list[i].normal = glm::normalize(glm::cross(triangle.vertices[1] - triangle.vertices[0],triangle.vertices[2] - triangle.vertices[0]));
	}
}

//calculate the map<point,std::vector<ModelTriangle>> which keys represent the point and the value is the list of ModelTriangle that containing that point
std::map<glm::vec3,std::vector<ModelTriangle>> calculateVertexTriangleMap(std::vector<ModelTriangle> triangles){

}

std::vector<float> interpolateSingleFloats(float start, float end, int steps) {
	std::vector<float> result;
	float step = (end - start) / float(steps-1);
	for (int i = 0; i < steps; i++) {
		result.push_back(start + i * step);
	}
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 start, glm::vec3 end, int steps) {
	std::vector<glm::vec3> result;
	glm::vec3 step = (end - start) / float(steps-1);
	for (int i = 0; i < steps; i++) {
		result.push_back(start + float(i) * step);
	}
	return result;
}


glm::vec3 readVertex(std::vector<std::string> line, std::vector<glm::vec3> vertices){
	return glm::vec3(std::stof(line[1]),std::stof(line[2]),std::stof(line[3]));
}

/**
	 * @brief For 4th week, I just remove the last character.
	 * @deprecated Only work for week 4's task
	 *
	 * @param line a line from io stream
	 * @param vertices vertices that previously store
	 */
ModelTriangle formTriangle(std::vector<std::string> line, std::vector<glm::vec3> vertices,Colour mtl){

	for (int i = 1; i < line.size(); i++){
		line[i].pop_back();
	}

	int i1, i2, i3;
	i1 = std::stoi(line[1])-1;
	i2 = std::stoi(line[2])-1;
	i3 = std::stoi(line[3])-1;


	return ModelTriangle(vertices[i1],vertices[i2],vertices[i3],mtl);
}

std::vector<ModelTriangle> readObj(std::string filepath,std::map<std::string,Colour> mtlmap,float scaling){
	std::string receiver;
	std::ifstream readFile(filepath);

	std::vector<glm::vec3> vertices;
	std::vector<ModelTriangle> triangles;

	Colour currentmtl;

	while(getline(readFile,receiver)){
		std::vector<std::string> elements = split(receiver,' ');
		if (elements[0] == "v"){
			vertices.push_back(readVertex(elements,vertices)*scaling);
		}else if(elements[0] == "f"){
			triangles.push_back(formTriangle(elements,vertices,currentmtl));
		}else if (elements[0] == "usemtl"){
			currentmtl = mtlmap[elements[1]];
		}
	}
	readFile.close();
	std::cout << "READ OBJ COMPLETE" << std::endl;
	return triangles;
}

void addmtl(std::vector<std::string> line, std::vector<std::string> elements, std::map<std::string,Colour> &map){
	map[line[1]] = Colour(line[1],std::stof(elements[1])*255,std::stof(elements[2])*255,std::stof(elements[3])*255);

}

std::map<std::string,Colour> readMtl(std::string filepath){
	std::string receiver;
	std::ifstream readFile(filepath);
	std::map<std::string,Colour> mtlmap;

	while(getline(readFile,receiver)){
		std::vector<std::string> line = split(receiver,' ');
		if(line[0] == "newmtl"){
			// std::cout << receiver << std::endl;
			getline(readFile,receiver);
			// std::cout << receiver << std::endl;
			std::vector<std::string> nextline = split(receiver,' ');
			addmtl(line,nextline,mtlmap);
		}
	}
	std::cout << "READ MTL COMPLETE" << std::endl;

	return mtlmap;
}


//Drawing Functions Start Here: 

//________________________________________ .................................____________________________________________
void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, DrawingWindow &window){
	float xdiff = to.x - from.x;
	float ydiff = to.y - from.y;
	float step = std::max(abs(xdiff),abs(ydiff));
	float xStepSize = xdiff/step;
	float yStepSize = ydiff/step;

	for (int i = 0; i < step; i++){
		float x = from.x+( xStepSize * i );
		float y = from.y + ( yStepSize * i );
		pbuffer[y][x] = fromColour(colour);
	}
}

void drawTriangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window){
	 drawLine(triangle.v0(),triangle.v1(),colour,window);
	 drawLine(triangle.v1(),triangle.v2(),colour,window);
	 drawLine(triangle.v2(),triangle.v0(),colour,window);
}

void drawFlatFilledTriangle(CanvasPoint bottom1, CanvasPoint bottom2, CanvasPoint top, Colour colour, DrawingWindow &window){
	float direction = top.y > bottom1.y ? 1 : -1;
	float slope1 = (top.x - bottom1.x)/(top.y-bottom1.y);
	float slope2 = (top.x - bottom2.x)/(top.y-bottom2.y);
	for (int i = 0; i < abs(top.y - bottom1.y); i++){
		CanvasPoint b1 = CanvasPoint(
		bottom1.x+ slope1*(i*direction),
		//x      + (      dx          /      dy        )*(    dy     )
		bottom1.y+(i*direction)
		);
		CanvasPoint b2 = CanvasPoint(
		bottom2.x+ slope2*(i*direction),
		bottom2.y+(i*direction)
		);
		drawLine(b1,b2,colour,window);
	}


}

// void drawFilledTriangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window){
// 	//Sort Points From top to Bottom
// 	CanvasPoint top,middle,bottom,split;
// 	std::array<CanvasPoint,3> points= triangle.vertices;
// 	std::sort(points.begin(),points.end(),[](CanvasPoint a, CanvasPoint b){return a.y<b.y;});

// 	top= points[0];
// 	middle= points[1];
// 	bottom=points[2];
// 	//Assertions
// 	assert(middle.y >= top.y);
//     assert(bottom.y >= middle.y);

// 	//Calculate Split Point
// 	split = CanvasPoint(
// 		top.x + (bottom.x - top.x)/(bottom.y-top.y)*(middle.y-top.y),
// 		middle.y
// 	);
// 	drawFlatFilledTriangle(middle,split,top,colour,window);
// 	drawFlatFilledTriangle(middle,split,bottom,colour,window);
	
// }

glm::mat3x3 transferFunction(CanvasTriangle triangle1, CanvasTriangle triangle2){
	// glm::mat3x3 mat1(triangle1.v0().x,triangle1.v0().y,1,triangle1.v1().x,triangle1.v1().y,1,triangle1.v2().x,triangle1.v2().y,1);
	// glm::mat3x3 mat2(triangle2.v0().x,triangle2.v0().y,1,triangle2.v1().x,triangle2.v1().y,1,triangle2.v2().x,triangle2.v2().y,1);
	glm::mat3x3 mat1(triangle1.v0().x,triangle1.v1().x,triangle1.v2().x,triangle1.v0().y,triangle1.v1().y,triangle1.v2().y,1,1,1);
	glm::mat3x3 mat2(triangle2.v0().x,triangle2.v1().x,triangle2.v2().x,triangle2.v0().y,triangle2.v1().y,triangle2.v2().y,1,1,1);
	glm::mat3x3 result = glm::inverse(mat1)*mat2;
	return result;
}

void drawTexturedLine(CanvasPoint from, CanvasPoint to, TextureMap texture, glm::mat3x3 transfer, DrawingWindow &window){
	float xdiff = to.x - from.x;
	float ydiff = to.y - from.y;
	float step = std::max(abs(xdiff),abs(ydiff));
	float xStepSize = xdiff/step;
	float yStepSize = ydiff/step;

	for (int i = 0; i < step; i++){
		float x = from.x+( xStepSize * i );
		float y = from.y + ( yStepSize * i );
		glm::vec3 p = glm::vec3(x,y,1)*transfer ;
		TexturePoint point(p.x,p.y);
		uint32_t c = texture.pixels[int(p.y)*texture.width+int(p.x)];
		window.setPixelColour(round(x),round(y),c);
	}
}

void drawFlatFilledTriangle(CanvasPoint bottom1, CanvasPoint bottom2, CanvasPoint top, TextureMap texture, glm::mat3x3 transfer, DrawingWindow &window){
	//detect drawing direction
	float directiony = top.y >= bottom1.y ? 1 : -1;

	float slope1 = (top.x - bottom1.x)/(top.y-bottom1.y);
	float slope2 = (top.x - bottom2.x)/(top.y-bottom2.y);
	int height = abs(top.y - bottom1.y);
	for (int i = 0; i < height; i++){
		float factor = i*directiony;
		CanvasPoint b1 = CanvasPoint(
		bottom1.x+ slope1*factor,
		//x      + (      dx          /      dy        )*(    dy     )
		bottom1.y+factor
		);
		CanvasPoint b2 = CanvasPoint(
			bottom2.x+slope2*factor,
			bottom2.y+factor
		);
		drawTexturedLine(b1,b2,texture,transfer,window);
		
	}
}

void drawFilledTriangle(CanvasTriangle triangle, DrawingWindow &window){
	CanvasTriangle texture = CanvasTriangle(CanvasPoint(triangle.v0().texturePoint.x,triangle.v0().texturePoint.y),
											CanvasPoint(triangle.v1().texturePoint.x,triangle.v1().texturePoint.y),
											CanvasPoint(triangle.v2().texturePoint.x,triangle.v2().texturePoint.y)
											);
	glm::mat3x3 transfer = transferFunction(triangle,texture);
	TextureMap tmap = TextureMap("./src/texture.ppm");

	

	CanvasPoint top,middle,bottom,split;
	std::array<CanvasPoint,3> points= triangle.vertices;
	std::sort(points.begin(),points.end(),[](CanvasPoint a, CanvasPoint b){return a.y<b.y;});

	top= points[0];
	middle= points[1];
	bottom=points[2];
	drawTriangle(triangle,Colour(255,255,255),window);
	//Assertions
	assert(middle.y >= top.y);
    assert(bottom.y >= middle.y);

	//Calculate Split Point
	float splitx = top.x + (bottom.x - top.x)/(bottom.y-top.y)*(middle.y-top.y);
	float splity = middle.y;
	glm::vec3 tPoint (glm::vec3(splitx,splity,1)*transfer);
	TexturePoint texturePoint = TexturePoint(tPoint.x,tPoint.y);
	split = CanvasPoint(
		splitx,
		splity
	);
	split.texturePoint = texturePoint;
	drawFlatFilledTriangle(middle,split,top,tmap,transfer,window);
	drawFlatFilledTriangle(middle,split,bottom,tmap,transfer,window);
	
}







///////_____________________________________________
bool equalPoint(glm::vec3 point1, glm::vec3 point2, float diff){
	if (glm::distance(point1,point2) < diff){
		return true;
	}
	return false;
}

RayTriangleIntersection getClosestIntersection(glm::vec3 campos, glm::vec3 direction, std::vector<ModelTriangle> triangles){
	RayTriangleIntersection result(glm::vec3(0,0,0), FLT_MAX, ModelTriangle(glm::vec3(0,0,0),glm::vec3(0,0,0),glm::vec3(0,0,0),Colour(0,0,0)),-1);
	for (int i = 0; i < triangles.size(); i++){
		glm::vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
		glm::vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
		glm::vec3 SPVector = campos - triangles[i].vertices[0];
		glm::mat3 DEMatrix(-direction, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		if (possibleSolution[1] >= 0.0 && possibleSolution[1] <= 1.0 &&
			possibleSolution[2] >= 0.0 && possibleSolution[2] <= 1.0 &&
			possibleSolution[1] + possibleSolution[2] <= 1.0 &&
			possibleSolution[0] >= 0.0){

			if (result.distanceFromCamera>possibleSolution.x){
			// std::cout << "get three " << std::endl;
			result.distanceFromCamera = possibleSolution.x;
			glm::vec3 intersectionPoint = triangles[i].vertices[0] + e0*possibleSolution.y + e1*possibleSolution.z;
			// glm::vec3 intersectionPiont2 = campos + direction*possibleSolution.x;
			// assert(intersectionPoint1 == intersectionPiont2);
			result.intersectionPoint = intersectionPoint;
			result.intersectedTriangle = triangles[i];
			result.triangleIndex = i;
			result.possibleSolution = possibleSolution;
		}

	}
			//std::cout << possibleSolution.x <<" " << possibleSolution.y <<" "  << possibleSolution.z <<" " <<std::endl;
			
		}

		
	return result;
}


glm::vec3 GetBarycentricCoord(CanvasPoint P1, CanvasPoint P2, CanvasPoint P3, CanvasPoint P)
{
	float u = ((P2.y - P3.y) * P.x + (P3.x - P2.x) * P.y + (P2.x * P3.y - P3.x * P2.y)) / ((P2.y - P3.y) * P1.x + (P3.x - P2.x) * P1.y + (P2.x * P3.y - P3.x * P2.y));
	float v = ((P1.y - P3.y) * P.x + (P3.x - P1.x) * P.y + (P1.x * P3.y - P3.x * P1.y)) / ((P1.y - P3.y) * P2.x + (P3.x - P1.x) * P2.y + (P1.x * P3.y - P3.x * P1.y));
	//float w = ((P1.y - P2.y) * P.x + (P2.x - P1.x) * P.y + (P1.x * P2.y - P2.x * P1.y)) / ((P1.y - P2.y) * P3.x + (P2.x - P1.x) * P3.y + (P1.x * P2.y - P2.x * P1.y));
	float w = 1 - u - v;
	// std::cout <<"SB" <<std::endl;
	return glm::vec3(u, v, w);
	// return glm::vec3(1, 1, 1);
}

float getdepth(float x, float y, CanvasTriangle triangle){
	auto v0 = triangle.v0();
	auto v1 = triangle.v1();
	auto v2 = triangle.v2();
	glm::vec3 barycord = GetBarycentricCoord(v0,v1,v2,CanvasPoint(round(x),round(y)));
	glm::vec3 target = glm::mat3x3(v0.x,v0.y,v0.depth,v1.x,v1.y,v1.depth,v2.x,v2.y,v2.depth)*barycord;

	return target.z;
}


void drawFilledTriangle(CanvasTriangle triangle, Colour colour,DrawingWindow &window){
	auto v0 = triangle.v0();
	auto v1 = triangle.v1();
	auto v2 = triangle.v2();
	int startx, starty, endx, endy;
	startx = std::min({v0.x,v1.x,v2.x});
	starty = std::min({v0.y,v1.y,v2.y});
	endx = std::max({v0.x,v1.x,v2.x});
	endy = std::max({v0.y,v1.y,v2.y});

	// std::cout << "BEFORE LOOP" << std::endl;

	for (int y = round(starty); y<round(endy); y++){

		for (int x = round(startx); x<round(endx); x++){
			if (!(x >= WIDTH || x < 0 || y >= HEIGHT || y < 0)){
				glm::vec3 barycord = GetBarycentricCoord(v0,v1,v2,CanvasPoint(round(x),round(y)));
				if (barycord.x >= 0 && barycord.x <= 1 && barycord.y >= 0 && barycord.y <= 1 && barycord.z >= 0 && barycord.z <= 1){
					// float depth = 1/(barycord*glm::mat3x3(v0.x,v0.y,v0.depth,v1.x,v1.y,v1.depth,v2.x,v2.y,v2.depth)).z;
					float depth = 1.0/(barycord*glm::mat3x3(v0.x,v1.x,v2.x,v0.y,v1.y,v2.y,v0.depth,v1.depth,v2.depth)).z;
					// float depth = 0.5;

				if (depth >= dbuffer[y][x]){
					// if (depth < dbuffer[0][0]){
					dbuffer[y][x] = depth;
					pbuffer[y][x] = fromColour(colour);
					//window.setPixelColour(x,y,fromColour(colour));
				}
			}
			}

		}
	}


}

glm::vec3 calculateDirection(float x, float y, glm::vec3 campos, glm::mat3 camrot){
	float focal = 2.0;
	float factor = 1.0/(HEIGHT*2.0/3.0);

	// float xres = (x-WIDTH/2)*(campos.z + focal)/(factor*focal);
	// float yres = (y+HEIGHT/2)*(campos.z + focal) /(factor*focal);
	float xres = (x-WIDTH/2.0) * factor;
	float yres = (y-HEIGHT/2.0) * factor;

	//assert(xres != 0 && yres != 0);
	glm::vec3 result(xres,yres,-focal);
	return camrot*result;

}


//Given a point in WORLD,camera position in WORLD, and the rotation ralative to WORLD and calculate the coreesponding render pixel in canvas
CanvasPoint calculateCanvasPos(glm::vec3 point,glm::vec3 campos,glm::mat3 camrot){
	//default focal length 2
	float focal = 2.0;

	//first move camera to world position, and then rotate.
	// glm::vec3 newPoint = glm::inverse(camrot) * glm::vec3(point.x - campos.x,point.y-campos.y,point.z-campos.z);
	glm::vec3 newPoint = glm::inverse(camrot) * glm::vec3(campos.x - point.x,campos.y-point.y,campos.z-point.z);

	float x = focal * ((newPoint.x)/(newPoint.z))*(-(HEIGHT*2.0/3.0)) + WIDTH/2;
	float y = focal * ((newPoint.y)/(newPoint.z))*(-(HEIGHT*2.0/3.0)) + HEIGHT/2;
	// printf("DBG:%f,%f POINT:%f,%f,%f\n",x,y,point.x,point.y,point.z);
	// return CanvasPoint(round(x),round(y),glm::length(camrot*glm::vec3(0,0,-1)*newPoint));
	return CanvasPoint(round(x),round(y),newPoint.z);


}




// void drawPoints(DrawingWindow &window,std::vector<ModelTriangle> triangles){
// 	for (int i = 0; i < triangles.size(); i++){
// 		for (int j = 0; j < 3; j++){
// 			CanvasPoint point = calculateCanvasPos(triangles[i].vertices[j],);
// 			window.setPixelColour(point.x,point.y,fromColour(Colour(rand()%255,rand()%255,rand()%255)));
// 		}
// 	}
// }





//Main Start Here

CanvasPoint generateRandPoint(){
	CanvasPoint point = CanvasPoint(rand()%WIDTH, rand()%HEIGHT);
	return point;
}

CanvasTriangle generateRandTriangle(){
	CanvasTriangle triangle(generateRandPoint(),generateRandPoint(),generateRandPoint());
	return triangle;
}

Colour generateRandColour(){
	return Colour(rand()%255,rand()%255,rand()%255);
}

void drawModelTriangle( ModelTriangle triangle, glm::vec3 campos, glm::mat3x3 rot, DrawingWindow &window){
	std::array<glm::vec3, 3UL> vertices = triangle.vertices;

	//rotation matrix

	// glm::mat3 rot(1,0,0,0,0.7,0.7,0,-0.7,0.7);
	// glm::mat3 rot(1,0,0,0,1,0,0,0,1);


	CanvasTriangle triangle2d = CanvasTriangle(
		calculateCanvasPos(vertices[0],campos,rot),
		calculateCanvasPos(vertices[1],campos,rot),
		calculateCanvasPos(vertices[2],campos,rot)
	);
	drawFilledTriangle(triangle2d,triangle.colour,window);
	//drawTriangle(triangle2d,Colour(255,255,255),window);
}

void drawWiredModelTriangle( ModelTriangle triangle, glm::vec3 campos, glm::mat3x3 rot, DrawingWindow &window){
	std::array<glm::vec3, 3UL> vertices = triangle.vertices;

	//rotation matrix

	// glm::mat3 rot(1,0,0,0,0.7,0.7,0,-0.7,0.7);
	// glm::mat3 rot(1,0,0,0,1,0,0,0,1);


	CanvasTriangle triangle2d = CanvasTriangle(
		calculateCanvasPos(vertices[0],campos,rot),
		calculateCanvasPos(vertices[1],campos,rot),
		calculateCanvasPos(vertices[2],campos,rot)
	);
	// drawFilledTriangle(triangle2d,triangle.colour,window);
	drawTriangle(triangle2d,Colour(255,255,255),window);
}



glm::mat3x3 lookat(glm::vec3 campos, glm::vec3 point=glm::vec3(0,0,0)){
	glm::vec3 fakeup(0,1,0);
	glm::vec3 forward = glm::normalize(-point + campos);
	glm::vec3 right = glm::normalize(glm::cross(fakeup,forward));
	glm::vec3 up = glm::normalize(glm::cross(right,forward));
	return glm::mat3(right,up,forward);
}


void drawRasterizedCameraView(DrawingWindow &window, glm::vec3 campos, glm::mat3x3 camrot, std::vector<ModelTriangle> list){

	for (int i =0; i < list.size(); i++){
		drawModelTriangle(list[i],campos,camrot,window);
		// std::cout << "Complete triangle " << i << std::endl;
	}

}


void drawWireFrameCameraView(DrawingWindow &window, glm::vec3 campos, glm::mat3x3 camrot, std::vector<ModelTriangle> list){
	for(ModelTriangle triangle : list){
		drawWiredModelTriangle(triangle,campos,camrot,window);
	}
}


void drawRaytracingCameraView(DrawingWindow &window, glm::vec3 campos, glm::mat3x3 camrot, std::vector<ModelTriangle> list, std::vector<glm::vec3> lightList){
	for( int y = 0 ; y < HEIGHT; y++){
		for ( int x = 0; x < WIDTH; x++){
			//Send Ray
			glm::vec3 raydirection = glm::normalize(calculateDirection(x,y,campos,camrot));
			RayTriangleIntersection camintersect =  getClosestIntersection(campos,raydirection,list);

			//If hit
			if (camintersect.triangleIndex != -1){
				Colour original = camintersect.intersectedTriangle.colour;
				float ambient = 0.0;
				float lightIntensity = 0.0;
				float diffuseIntensity = 0.0;
				float specularIntensity = 0.0;
				float visibility = 1.0;
				for(int i = 0; i < lightList.size(); i++){

					glm::vec3 lights = lightList[i];
					RayTriangleIntersection lightintersect = getClosestIntersection(lights,camintersect.intersectionPoint-lights,list);
					if(equalPoint(lightintersect.intersectionPoint, camintersect.intersectionPoint,0.01)){
						ModelTriangle triangle = lightintersect.intersectedTriangle;
						glm::vec3 vectorToLight = glm::normalize(lights - lightintersect.intersectionPoint);
						glm::vec3 vectorFromLight = glm::normalize( lightintersect.intersectionPoint- lights);
						glm::vec3 reflection = vectorFromLight- 2.0f * triangle.normal*(glm::dot(vectorFromLight,triangle.normal));
						
						
						// float normalIntensity = glm::length(triangle.normal*vectorToLight);
						lightIntensity += (1/(4*PI*pow(lightintersect.distanceFromCamera,2.0)));
						diffuseIntensity += glm::dot(triangle.normal,vectorToLight);
						specularIntensity += pow(std::max(glm::dot(reflection,glm::normalize(campos-lightintersect.intersectionPoint)),0.0f),256);
						// std::cout << specularIntensity << std::endl;
						
						visibility = glm::dot((campos-lightintersect.intersectionPoint),triangle.normal)+visibility > 0 ? 1.0 : 0.0;
						
						// std::cout << intensity << std::endl;
						// assert(triangle.normal != glm::vec3(0,0,0));
					
					}	
				}
				lightIntensity /= float(lightList.size());
				diffuseIntensity /= float(lightList.size());
				specularIntensity /= float(lightList.size());




				glm::vec3 colour = glm::vec3(original.red,original.green,original.blue);
					
				glm::vec3 resultColor = colour*(lightIntensity+diffuseIntensity+specularIntensity+ambient)*visibility;
				// glm::vec3 resultColor = colour*(specularIntensity);



				// float r = std::min(float(255),original.red*(intensity+amb+normalIntensity+specularIntensity)*visibility);
				// float g = std::min(float(255),original.green*(intensity+amb+normalIntensity+specularIntensity)*visibility);
				// float b = std::min(float(255),original.blue*(intensity+amb+normalIntensity+specularIntensity)*visibility);


				Colour pixel = Colour(std::min(resultColor.x,255.0f),std::min(resultColor.y,255.0f),std::min(resultColor.z,255.0f));

				// pbuffer[y][x] = fromColour(camintersect.intersectedTriangle.colour);
				pbuffer[y][x] = fromColour(pixel);
			}
			
		}
	}

}


glm::mat3x3 xzrot(float angle){
	return glm::mat3x3(cos(angle),0,sin(angle),0,1,0,-sin(angle),0,cos(angle));
}


void testReadObj(){
	printf("print result \n");
	auto mtlmap = readMtl("./src/cornell-box.mtl");
	auto list = readObj("./src/cornell-box.obj",mtlmap,1);
	for(int i = 0; i < list.size(); i++){
		std::cout<<list[i]<<std::endl;
		std::cout<<list[i].colour << std::endl;
	}
	printf("print all %lu result\n",list.size());

}




void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) campos.x-=0.1;
		else if (event.key.keysym.sym == SDLK_RIGHT) campos.x+=0.1;
		else if (event.key.keysym.sym == SDLK_UP) campos.y-=0.1;
		else if (event.key.keysym.sym == SDLK_DOWN)campos.y+=0.1;
		else if (event.key.keysym.sym == SDLK_RIGHTBRACKET)campos.z+=0.1;
		else if (event.key.keysym.sym == SDLK_LEFTBRACKET)campos.z-=0.1;
		else if (event.key.keysym.sym == SDLK_p)std::cout << glm::to_string(campos) << std::endl;
		else if (event.key.keysym.sym == SDLK_o){
			drawTriangle(generateRandTriangle(),generateRandColour(),window);
		}
		else if (event.key.keysym.sym == SDLK_i){
			CanvasTriangle triangle(
				CanvasPoint(160,10),
				CanvasPoint(300,230),
				CanvasPoint(10,150)
			);
			triangle.v0().texturePoint = TexturePoint(195,5);
			triangle.v1().texturePoint = TexturePoint(395,380);
			triangle.v2().texturePoint = TexturePoint(65,330);



			drawFilledTriangle(triangle,window);
		}
		else if (event.key.keysym.sym == SDLK_r){
			state = 0x02;
		}
		else if (event.key.keysym.sym == SDLK_e){
			state = 0x01;
		}
		else if (event.key.keysym.sym == SDLK_w){
			state = 0x03;
		}



		else if (event.type == SDL_MOUSEBUTTONDOWN) {
			window.savePPM("output.ppm");
			window.saveBMP("output.bmp");
		}
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	// drawFilledTriangle(generateRandTriangle(),generateRandColour(),window);
	std::map<std::string,Colour> mtl = readMtl("src/cornell-box.mtl");
	// std::vector<ModelTriangle> list = readObj("src/cornell-box.obj",mtl,0.35);
	std::vector<ModelTriangle> list = readObj("src/ball.obj",mtl,0.5);

	calculateNormal(list);
	state = 0x02;
	// for (size_t i = 0; i < list.size(); i++){

	// std::cout << list[i] << std::endl;
	// }

	campos = glm::vec3(0.000000, 0.800000, 2.200002);
	// campos = glm::vec3(0, 0, 4);

	// float angle = 0.01;
	// glm::mat3x3 rot(cos(angle),0,sin(angle),0,1,0,-sin(angle),0,cos(angle));
	while (true) {
		// angle += 0.001;
		glm::mat3 rot(1,0,0,0,-1,0,0,0,1);
		// glm::mat3x3 rot(cos(angle),0,sin(angle),0,1,0,-sin(angle),0,cos(angle));
		// glm::mat3x3 rotpos (cos(angle),0,sin(angle),0,1,0,-sin(angle),0,cos(angle));
		// campos = campos*rotpos;
		
		// glm::mat3x3 rot = lookat(campos);
		
		

		// std::cout << glm::to_string(rot)<< std::endl; 
		// std::cout << glm::to_string(rot) << std::endl;
		

		window.clearPixels();
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		
		if (state == 0x01){
			drawRasterizedCameraView(window,campos,rot,list);
		}
		if (state == 0x02){
			std::vector<glm::vec3> lights = {glm::vec3(-0.4,0.7,0.9)};
			drawRaytracingCameraView(window,campos,rot,list ,lights);
		}
		if (state == 0x03){
			drawWireFrameCameraView(window,campos,rot,list);
		}
		
		

		// std::cout << "campos: " << campos.x << campos.y << campos.z << std::endl;


		//END OF THE LOOP
		for (int y = 0; y < HEIGHT; y++){
			for (int x = 0; x < WIDTH ; x++){
				window.setPixelColour(x,y,pbuffer[y][x]);
				// std::cout<< pbuffer[y][x] << std::endl;
			}
		}
		dbuffer = makedBuffer();
		pbuffer = makepBuffer();
		window.renderFrame();
	}
}
