#include <iostream>
#include <GLFW/glfw3.h>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <random>
#include <algorithm>

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
    vec2 vel = vec2(0.0f);
	vec2 force = vec2(0.0f);

    float mass = 2;
    float charge = 1.0f;

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
	Atom(vec2(200.0f, 0.0f)),

    Atom(vec2(200.0f, 200.0f)),
    Atom(vec2(200.0f, -200.0f)),

    Atom(vec2(-200.0f, 200.0f)),
    Atom(vec2(-200.0f, -200.0f)),

    Atom(vec2(0.0f, 200.0f)),
    Atom(vec2(0.0f, -200.0f)),
};


// ================= Orbital Sampling ================= //
const double a0 = 1.0; // Bohr radius

// Quantum numbers of the orbital being sampled
int n = 2, l = 1, m = 0;

mt19937 gen(random_device{}());

// Associated Laguerre polynomial L_{n-l-1}^{2l+1}(rho)
double laguerre(int n, int l, double rho) {
    int k = n - l - 1;
    int alpha = 2 * l + 1;

    double L = 1.0, Lm1 = 1.0 + alpha - rho;
    if (k == 1) L = Lm1;
    else if (k > 1) {
        double Lm2 = 1.0;
        for (int j = 2; j <= k; ++j) {
            L = ((2 * j - 1 + alpha - rho) * Lm1 - (j - 1 + alpha) * Lm2) / j;
            Lm2 = Lm1;
            Lm1 = L;
        }
    }
    return L;
}

// Associated Legendre polynomial P_l^m(x)
double legendre(int l, int m, double x) {
    double Pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int j = 1; j <= m; ++j) {
            Pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m) return Pmm;

    double Pm1m = x * (2 * m + 1) * Pmm;
    if (l == m + 1) return Pm1m;

    for (int ll = m + 2; ll <= l; ++ll) {
        double Pll = ((2 * ll - 1) * x * Pm1m - (ll + m - 1) * Pmm) / (ll - m);
        Pmm = Pm1m;
        Pm1m = Pll;
    }
    return Pm1m;
}

// Radial distance, sampled from r^2 * R(r)^2 using a CDF
double sampleR(int n, int l) {
    const int N = 4096;
    const double rMax = 10.0 * n * n * a0;

    // The table is rebuilt whenever n or l changes
    static vector<double> cdf;
    static int builtN = -1, builtL = -1;

    if (builtN != n || builtL != l) {
        cdf.resize(N);
        double dr = rMax / (N - 1);
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double r = i * dr;
            double rho = 2.0 * r / (n * a0);

            double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));
            double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * laguerre(n, l, rho);

            sum += r * r * R * R;
            cdf[i] = sum;
        }

        for (double& v : cdf) v /= sum;
        builtN = n;
        builtL = l;
    }

    uniform_real_distribution<double> dis(0.0, 1.0);
    int idx = lower_bound(cdf.begin(), cdf.end(), dis(gen)) - cdf.begin();
    return idx * (rMax / (N - 1));
}

// Polar angle, sampled from sin(theta) * P_l^m(cos(theta))^2 using a CDF
double sampleTheta(int l, int m) {
    const int N = 2048;

    static vector<double> cdf;
    static int builtL = -1, builtM = 0;

    if (builtL != l || builtM != m) {
        cdf.resize(N);
        double dtheta = M_PI / (N - 1);
        double sum = 0.0;

        for (int i = 0; i < N; ++i) {
            double theta = i * dtheta;
            double Plm = legendre(l, abs(m), cos(theta));

            sum += sin(theta) * Plm * Plm;
            cdf[i] = sum;
        }

        for (double& v : cdf) v /= sum;
        builtL = l;
        builtM = m;
    }

    uniform_real_distribution<double> dis(0.0, 1.0);
    int idx = lower_bound(cdf.begin(), cdf.end(), dis(gen)) - cdf.begin();
    return idx * (M_PI / (N - 1));
}

// Azimuth is uniform since |e^(i m phi)|^2 = 1
double samplePhi() {
    uniform_real_distribution<double> dis(0.0, 1.0);
    return 2.0 * M_PI * dis(gen);
}

// y is up, theta is measured from the y axis
vec3 sphericalToCartesian(double r, double theta, double phi) {
    return vec3(r * sin(theta) * cos(phi),
                r * cos(theta),
                r * sin(theta) * sin(phi));
}

vector<vec3> generateOrbitalPoints(int count) {
    vector<vec3> points;
    points.reserve(count);
    for (int i = 0; i < count; ++i)
        points.push_back(sphericalToCartesian(sampleR(n, l), sampleTheta(l, m), samplePhi()));
    return points;
}


// ================= Orbital View ================= //

// Press O to switch between the atoms and the sampled orbital
bool orbitalView = false;

// Orbits around the nucleus: drag to rotate, scroll to zoom
struct Camera {
    float radius = 50.0f;
    float azimuth = 0.6f;
    float elevation = 1.2f; // measured from the y axis
    bool dragging = false;
    double lastX = 0.0, lastY = 0.0;

    vec3 position() const {
        return vec3(radius * sin(elevation) * cos(azimuth),
                    radius * cos(elevation),
                    radius * sin(elevation) * sin(azimuth));
    }

    void processMouseMove(double x, double y) {
        if (dragging) {
            azimuth += float(x - lastX) * 0.01f;
            elevation = std::clamp(elevation - float(y - lastY) * 0.01f, 0.01f, float(M_PI) - 0.01f);
        }
        lastX = x;
        lastY = y;
    }
};
Camera camera;

// black -> purple -> red -> orange -> yellow -> white
vec3 heatmapFire(float value) {
    const vec3 colors[6] = {
        vec3(0.0f, 0.0f, 0.0f),
        vec3(0.5f, 0.0f, 0.99f),
        vec3(0.8f, 0.0f, 0.0f),
        vec3(1.0f, 0.5f, 0.0f),
        vec3(1.0f, 1.0f, 0.0f),
        vec3(1.0f, 1.0f, 1.0f),
    };

    float scaled = std::clamp(value, 0.0f, 1.0f) * 5.0f;
    int i = std::min((int)scaled, 4);
    return mix(colors[i], colors[i + 1], scaled - i);
}

// |R(r)|^2 * |P_l^m(cos(theta))|^2 at a point, the same weight the points were sampled from
double orbitalDensity(const vec3& p) {
    double r = length(p);
    if (r < 1e-6) return 0.0;

    double rho = 2.0 * r / (n * a0);
    double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));
    double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * laguerre(n, l, rho);
    double P = legendre(l, abs(m), p.y / r);
    return R * R * P * P;
}

// Colors each point by its density, scaled so the denser half of the cloud shows up as warm colors
vector<vec3> colorOrbitalPoints(const vector<vec3>& points) {
    vector<double> density;
    density.reserve(points.size());
    for (const vec3& p : points)
        density.push_back(orbitalDensity(p));

    vector<double> sorted = density;
    sort(sorted.begin(), sorted.end());
    double ref = std::max(sorted[(size_t)(sorted.size() * 0.6)], 1e-30);

    vector<vec3> colors;
    colors.reserve(points.size());
    for (double d : density)
        colors.push_back(heatmapFire((float)pow(d / ref, 0.6)));
    return colors;
}

// Flat grid on the y = 0 plane
void drawGrid(float half) {
    glColor3f(0.6f, 0.6f, 0.6f);
    glBegin(GL_LINES);
    for (int i = -1; i <= 1; i++) {
        float c = i * half;
        glVertex3f(-half, 0.0f, c); glVertex3f(half, 0.0f, c);
        glVertex3f(c, 0.0f, -half); glVertex3f(c, 0.0f, half);
    }
    glEnd();
}

// Draws the points in 3D from the orbit camera, looking at the nucleus at the origin
void drawOrbitalPoints(const vector<vec3>& points, const vector<vec3>& colors, int width, int height) {
    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    mat4 projection = perspective(radians(45.0f), width / (float)height, 0.1f, 2000.0f);
    mat4 view = lookAt(camera.position(), vec3(0.0f), vec3(0.0f, 1.0f, 0.0f));

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(value_ptr(projection));
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(value_ptr(view));

    drawGrid(250.0f);

    glEnable(GL_POINT_SMOOTH);
    glPointSize(2.0f * height / 600.0f);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < points.size(); i++) {
        glColor3f(colors[i].r, colors[i].g, colors[i].b);
        glVertex3f(points[i].x, points[i].y, points[i].z);
    }
    glEnd();

    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_DEPTH_TEST);
}


int main()
{
    Engine engine;

    glfwSetKeyCallback(engine.window, [](GLFWwindow*, int key, int, int action, int) {
        if (key == GLFW_KEY_O && action == GLFW_PRESS)
            orbitalView = !orbitalView;
    });
    glfwSetMouseButtonCallback(engine.window, [](GLFWwindow* win, int button, int action, int) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            camera.dragging = (action == GLFW_PRESS);
            glfwGetCursorPos(win, &camera.lastX, &camera.lastY);
        }
    });
    glfwSetCursorPosCallback(engine.window, [](GLFWwindow*, double x, double y) {
        camera.processMouseMove(x, y);
    });
    glfwSetScrollCallback(engine.window, [](GLFWwindow*, double, double yoffset) {
        camera.radius = std::max(1.0f, camera.radius - float(yoffset) * 0.1f * camera.radius);
    });

    float k = 1000.0f;
	float epsilon = 50.0f;

    // Sample the orbital once, it is drawn when O is pressed
    vector<vec3> orbitalPoints = generateOrbitalPoints(50000);
    float meanR = 0.0f;
    for (vec3& p : orbitalPoints)
        meanR += length(p);
    cout << "Sampled " << orbitalPoints.size() << " points for n=" << n << " l=" << l << " m=" << m
         << ", mean r = " << meanR / orbitalPoints.size() << endl;

    vector<vec3> orbitalColors = colorOrbitalPoints(orbitalPoints);
    camera.radius = 4.0f * meanR / orbitalPoints.size();

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(engine.window))
    {
		//Projection to Orthographic -1 to 1 --> 0 to 800, 0 to 600
        engine.run();

        if (orbitalView) {
            int w, h;
            glfwGetFramebufferSize(engine.window, &w, &h);
            drawOrbitalPoints(orbitalPoints, orbitalColors, w, h);
            glfwSwapBuffers(engine.window);
            glfwPollEvents();
            continue;
        }

		// Calculate forces
        for (Atom& a : atoms)
            a.force = vec2(0.0f);

        for (size_t i = 0; i < atoms.size(); i++) {
            for (size_t j = i + 1; j < atoms.size(); j++) {
                vec2 dir = atoms[j].pos - atoms[i].pos;
				float dist = length(dir);

				if (dist < 1.0f) dist = 1.0f; // Prevent division by zero
				vec2 unitdir = normalize(dir);

                // Coulumb Calc
				float forcemag = k * (atoms[i].charge * atoms[j].charge) / (dist * dist + epsilon);
                vec2 force = forcemag * unitdir;

                atoms[i].force += force;
                atoms[j].force += force;
            }
        }

        float dt = 0.016f;

        for (Atom& a : atoms) {
            vec2 acc = a.force / a.mass;
            a.vel += acc * dt;
            a.pos += a.vel * dt;
        }

		// clear the screen

        glClear(GL_COLOR_BUFFER_BIT);
        /* Render here */

        for (Atom& a : atoms) {
            for (Particle& p : a.particles) {
                if (p.charge == 1) {
                    p.pos = a.pos; // Update the center of the particle to the current position of the atom     
                }; 
                if (p.charge == -1) {
                    p.center = a.pos;
                    p.update();
                }

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