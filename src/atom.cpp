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
	vec2 center;
	vec2 velocity = vec2(0.0f);
    vec2 force = vec2(0.0f);

    int charge;
    float mass = 0.1f;

    float angle = 0.0f;

    Particle(vec2 pos, int charge, vec2 center) : pos(pos), charge(charge), center(center) {}

    void draw(int segments = 50) {
        float r;

        if (charge != -1) {
            glLineWidth(0.01f);
            glBegin(GL_LINE_LOOP);
            glColor3f(0.4f, 0.4f, 0.4f); // grey for the orbit
            for (int i = 0; i <= segments; i++) {
                float angle = 2.0f * M_PI * i / segments;
                float x = cos(angle) * orbitDistance;
                float y = sin(angle) * orbitDistance;
                glVertex2f(x + pos.x, y + pos.y);
            }
            glEnd();
        }


        if (charge == -1) {
            r = 5;
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
        angle += 0.0025f; // Increment the angle for rotation
		pos = center + vec2(cos(angle) * orbitDistance,
                            sin(angle) * orbitDistance); // Update position to create circular motion)
    }
};

struct Atom {
    vec2 pos;
	vector<Particle> particles;

    Atom(vec2 pos) : pos(pos) {
		particles.emplace_back(pos, 1, pos); // Proton at the center
        particles.emplace_back(
			pos + vec2(orbitDistance,0), -1, pos); // Electron at the initial position on the orbit
    }
};

vector<Atom> atoms{
	Atom(vec2(0.0f)),
	Atom(vec2(-200.0f, 0.0f)),
};


int main()
{
    Engine engine;
    float k = 500.0f;

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(engine.window))
    {
		//Projection to Orthographic -1 to 1 --> 0 to 800, 0 to 600
        engine.run();

		// clear the screen

        glClear(GL_COLOR_BUFFER_BIT);
        /* Render here */

        for (Atom& a : atoms) {
            for (Particle& p : a.particles) {
				p.force = vec2(0.0f);
                if (p.charge == -1)
                    p.update();
            }
            for (Particle& p : a.particles) {
                p.draw();
            }
        }

        

        
        /* Poll for and process events */
        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    glfwDestroyWindow(engine.window);
    glfwTerminate();

}