#include <iostream>
#include <GLFW/glfw3.h>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace glm;

// Variables
float orbitDistance = 50.0f;


struct Engine {

    GLFWwindow* window;
	int WIDTH = 800, HEIGHT=600;

    void run() {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-WIDTH / 2.0, WIDTH / 2.0,
            -HEIGHT / 2.0, HEIGHT / 2.0,
            -1, 1);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }

    Engine(){
        /* Initialize the library */
        if (!glfwInit()) {
            cerr<< "failed to initialize GLFW" << endl;
             exit(EXIT_FAILURE);
        }
        /* Create a windowed mode window and its OpenGL context */
        window = glfwCreateWindow(WIDTH, HEIGHT, "Atom Simulation", NULL, NULL);
        if (!window)
        {
			cerr << "failed to create GLFW window" << endl;
            glfwTerminate();
        }
        /* Make the window's context current */
        glfwMakeContextCurrent(window);
	}

};

struct Particle {
    vec2 pos;
    int charge;
    float angle = 0.0f;

    Particle(vec2 pos, int charge) : pos(pos), charge(charge) {}

    void draw(int segments = 50) {
        float r;

        if (charge == -1) {
            r = 2;
            glColor3f(0.0f, 0.0f, 1.0f); // blue for electrons
        }
        else if (charge == 1) {
            r = 10;
            glColor3f(1.0f, 0.0f, 0.0f); // red for protons
        }
        else {
            r = 10;
            glColor3f(1.0f, 1.0f, 1.0f); // white for neutrons
        }

        // draw the particle
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(pos.x, pos.y); // center of the circle
        for (int i = 0; i <= segments; i++) {
            float angle = 2.0f * M_PI * i / segments;
            float x = cos(angle) * r;
            float y = sin(angle) * r;
            glVertex2f(x + pos.x, y + pos.y);
        }
        glEnd();
    };
    void update() {
        angle += 0.001f; // Increment the angle for rotation
		pos = vec2(cos(angle) * orbitDistance, sin(angle) * orbitDistance); // Update position to create circular motion)
    }
};



vector<Particle> particles = {
    Particle(vec2(0.0, 0.0), 1),
    Particle(vec2(-50.0, 0.0), -1)
};


int main()
{
    Engine engine;

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(engine.window))
    {
		//Projection to Orthographic -1 to 1 --> 0 to 800, 0 to 600
        engine.run();

		// clear the screen

        glClear(GL_COLOR_BUFFER_BIT);
        /* Render here */
        for (Particle& p : particles) {
            if (p.charge == -1)
                p.update();
        }

        for (Particle& p : particles) {
            p.draw();
        }
        /* Poll for and process events */
        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    glfwDestroyWindow(engine.window);
    glfwTerminate();

}