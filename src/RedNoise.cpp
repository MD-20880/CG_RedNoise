#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <DrawingWindow.h>
#include <TextureMap.h>
#include <Colour.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>	

#define WIDTH 320
#define HEIGHT 240






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

uint32_t fromColour(Colour colour){

	uint32_t colour_in_uint32 = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
	return colour_in_uint32;
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
		window.setPixelColour(round(x),round(y),fromColour(colour));
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

void drawFilledTriangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window){
	//Sort Points From top to Bottom
	CanvasPoint top,middle,bottom,split;
	std::array<CanvasPoint,3> points= triangle.vertices;
	std::sort(points.begin(),points.end(),[](CanvasPoint a, CanvasPoint b){return a.y<b.y;});

	top= points[0];
	middle= points[1];
	bottom=points[2];
	//Assertions
	assert(middle.y >= top.y);
    assert(bottom.y >= middle.y);

	//Calculate Split Point
	split = CanvasPoint(
		top.x + (bottom.x - top.x)/(bottom.y-top.y)*(middle.y-top.y),
		middle.y
	);
	drawFlatFilledTriangle(middle,split,top,colour,window);
	drawFlatFilledTriangle(middle,split,bottom,colour,window);
	
}

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
		glm::vec3 p = transfer * glm::vec3(x,y,1);
		TexturePoint point(p.x,p.y);
		uint32_t c = texture.pixels[int(p.y)*texture.width-1+int(p.x)];
		std::cout << "color at position " << p.x << " " << p.y << " is: "<< c << std::endl;
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



void drawFilledTriangle(CanvasTriangle triangle, CanvasTriangle texture, DrawingWindow &window){
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





void drawgray(DrawingWindow &window) {
	window.clearPixels();
	std::vector<float> grayscale = interpolateSingleFloats(255,0,WIDTH);
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = grayscale[x];
			float green = grayscale[x];
			float blue = grayscale[x];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}




void drawcolor(DrawingWindow &window) {
	window.clearPixels();
	glm::vec3 topleft(255,0,0);
	glm::vec3 topright(0,255,0);
	glm::vec3 bottomleft(0,0,255);
	glm::vec3 bottomright(255,255,0);

	std::vector<glm::vec3> firstCol = interpolateThreeElementValues(topleft,bottomleft,HEIGHT);
	std::vector<glm::vec3> lastCol = interpolateThreeElementValues(topright,bottomright,HEIGHT);

	for (int y = 0; y < HEIGHT; y++){
		std::vector<glm::vec3> row = interpolateThreeElementValues(firstCol[y],lastCol[y],WIDTH);
		for(int x = 0; x < WIDTH; x++){
			float red = row[x].x;
			float green = row[x].y;
			float blue = row[x].z;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
 

}






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




void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u){
			drawFilledTriangle(generateRandTriangle(),generateRandColour(),window);
		}
		else if (event.key.keysym.sym == SDLK_o){
			drawTriangle(generateRandTriangle(),generateRandColour(),window);
		}
		else if (event.key.keysym.sym == SDLK_i){
			CanvasTriangle triangle(
				CanvasPoint(160,10),
				CanvasPoint(300,230),
				CanvasPoint(10,150)
			);
			CanvasTriangle texture(
				CanvasPoint(195,5),
				CanvasPoint(395,380),
				CanvasPoint(65,330)
			);

			glm::mat3x3 matrix =  transferFunction(triangle,texture);
			glm::vec3 testvec(160,10,1);
			glm::vec3 result = testvec*matrix;
			std::cout << result.x << " " << result.y<< std::endl;

			drawFilledTriangle(triangle,texture,window);
			}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	} 
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// drawcolor(window);
		window.renderFrame();
	}
}
