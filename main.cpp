#include <fstream>
#include <iostream>
#include <SFML/Graphics.hpp>
#include <sstream>
#include <filesystem>
#include "include/Ensemble.h"
#include "include/QTEnsemble.h"
#include "include/VWParticle.h"
#include "include/Constants.h"



// Function declarations
void getEnsembleParameters(int& rate, int& iterations, float& dt, float& temperature, float& mass, bool& ideal, float& epsilon, float& sigma);
bool readSettingsFromCSV(const std::string& filename, int& rate, int& iterations, float& dt, float& temperature, float& mass, bool& ideal, float& epsilon, float& sigma);
double sampleMaxwellian(float T, float m);
std::vector<std::unique_ptr<Particle>> getRateParticles(Quad bounds, float dt, int rate, bool ideal, float T, float m, float epsilon = 0, float sigma = 0);
void createWalls(Wall walls[]);
void saveData(const double lT[], const double lP[], const double rT[], const double rP[], const double eT[], int iters, float dt);

int main() {
    // Ensemble parameters
    int rate, iterations;
    float dt, temperature, mass, epsilon, sigma;
    bool ideal;

    std::filesystem::path currentPath = std::filesystem::current_path();
    std::cout << currentPath << '\n';
    const bool fromCSV = readSettingsFromCSV(currentPath.string() + "/settings.csv", rate, iterations, dt, temperature, mass, ideal, epsilon, sigma);

    if (!fromCSV)
    {
        std::cout << "No settings file found, please enter the ensemble parameters manually." << std::endl;
        // Get ensemble parameters from the user
        getEnsembleParameters(rate, iterations, dt, temperature, mass, ideal, epsilon, sigma);
    }


    // Create walls for the simulation
    Wall ensembleBounds[14];
    createWalls(ensembleBounds);

    // Generate particles
    std::vector<std::unique_ptr<Particle>> particles;

    // Create the ensemble
    Quad QTBounds = Quad(WIDTH/2, HEIGHT/2, WIDTH, HEIGHT);
    QTEnsemble ensemble(std::move(particles), ensembleBounds, 14, QTBounds);

    // Create the SFML window
    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Particle Simulation");

    int i = 0;
    double leftChamberT[iterations], rightChamberT[iterations], ensembleT[iterations];
    double leftChamperP[iterations], rightChamberP[iterations];

    Quad spawnArea = Quad(MARGIN, MARGIN, (WIDTH - 2*MARGIN - THROTTLE_LENGTH)/2, HEIGHT - 2*MARGIN);
    Quad leftChamber = Quad(0, 0, WIDTH / 2, HEIGHT);
    Quad rightChamber = Quad(WIDTH / 2, 0, WIDTH / 2, HEIGHT);

    const int drawFrequency = 10;  // Draw every 10 iterations

    while (window.isOpen() && i < iterations) {
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Simulation logic
        ensemble.addParticles(getRateParticles(spawnArea, dt, rate, ideal, temperature, mass, epsilon, sigma));

        if (!ensemble.isEmpty()) {
            ensemble.cullNotInRegion(QTBounds);
            ensemble.iterateParticles(dt);

            ensembleT[i] = ensemble.getTemperature();
            leftChamberT[i] = ensemble.getTemperatureInRegion(leftChamber);
            rightChamberT[i] = ensemble.getTemperatureInRegion(rightChamber);

            leftChamperP[i] = ensemble.getPressureInRegion(leftChamber);
            rightChamberP[i] = ensemble.getPressureInRegion(rightChamber);
        }

        if (i % drawFrequency == 0) {
            // Clear, draw, and display only at specific intervals
            window.clear();
            ensemble.draw(window);
            window.display();
        }

        // Print progress bar every 10%
        if (i % (iterations / 10) == 0) {
            int progress = (i * 100 / iterations) / 10;  // Progress in 10% increments
            std::string bar = "[" + std::string(progress, '#') + std::string(10 - progress, ' ') + "]";
            std::cout << "\rProgress: " << bar << " " << progress * 10 << "%" << std::flush;
        }

        i++;
    }

    saveData(leftChamberT, leftChamperP, rightChamberT, rightChamberP, ensembleT, iterations, dt);

    return 0;
}

// Function to get ensemble parameters from the user
void getEnsembleParameters(int& rate, int& iterations, float& dt, float& temperature, float& mass, bool& ideal, float& epsilon, float& sigma) {
    std::cout << "Enter the rate of particle entry: ";
    std::cin >> rate;
    std::cout << "Enter the number of iterations [~10000]: ";
    std::cin >> iterations;
    std::cout << "Enter the time step [~0.01]: ";
    std::cin >> dt;
    std::cout << "Enter the temperature [K]: ";
    std::cin >> temperature;
    std::cout << "Enter the molar mass of particles [g mol^-1]: ";
    float molarMass;
    std::cin >> molarMass;
    mass = molarMass / 6.022e23 * 1e-3; // Convert molar mass to kg
    std::cout << "Is the ensemble ideal? (1 for yes, 0 for no): ";
    std::cin >> ideal;
    if (!ideal) {
        std::cout << "Enter the epsilon / k_b value: ";
        float reduced_epsilon;
        std::cin >> reduced_epsilon;
        epsilon = reduced_epsilon * K_b;
        std::cout << "Enter the sigma value [A]: ";
        float reduced_sigma;
        std::cin >> reduced_sigma;
        sigma = reduced_sigma * 1e-10; // Convert sigma to m
    }
}

// Function to read settings from a CSV file
bool readSettingsFromCSV(const std::string& filename, int& rate, int& iterations, float& dt, float& temperature, float& mass, bool& ideal, float& epsilon, float& sigma) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    std::cout << "== FOUND SETTINGS ==" << std::endl;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string key, value;
        if (std::getline(ss, key, ',') && std::getline(ss, value)) {
            std::cout << key << " = " << value << std::endl;
            if (key == "rate") rate = std::stoi(value);
            else if (key == "iterations") iterations = std::stoi(value);
            else if (key == "dt") dt = std::stof(value);
            else if (key == "temperature") temperature = std::stof(value);
            else if (key == "mass") mass = std::stof(value);
            else if (key == "ideal") ideal = std::stoi(value);
            else if (key == "epsilon") epsilon = std::stof(value);
            else if (key == "sigma") sigma = std::stof(value);
        }
    }
    file.close();
    return true;
}

// Function to sample a random speed from a 2D Maxwell-Boltzmann distribution
double sampleMaxwellian(const float T, const float m) {
    const double kB = 1.38064852e-23; // Boltzmann constant in J/K

    // Generate two independent random numbers uniformly distributed between 0 and 1
    double u1 = static_cast<double>(rand()) / RAND_MAX;
    double u2 = static_cast<double>(rand()) / RAND_MAX;

    // Use the Box-Muller transform to generate a normally distributed random number
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

    // Calculate the most probable speed (v_mp) for the 2D Maxwell-Boltzmann distribution
    double v_mp = sqrt(2.0 * kB * T / m);

    // Scale the normally distributed random number by the most probable speed
    double v = v_mp * z0;

    // Return the absolute value of the speed to ensure it's non-negative
    return fabs(v);
}

// Function to generate particles dynamically (Ideal or VW)
std::vector<std::unique_ptr<Particle>> generateParticles(int ensembleSize, bool ideal, float T, float m, float epsilon, float sigma) {
    std::vector<std::unique_ptr<Particle>> particles;
    particles.reserve(ensembleSize);
    srand(time(nullptr));
    float minDistance = 20;

    for (int i = 0; i < ensembleSize; i++) {
        bool validPosition = false;
        double x, y, vx, vy;

        while (!validPosition) {
            x = MARGIN + 10 + (std::rand() % (WIDTH / 3 - 20 - MARGIN));
            y = MARGIN + 10 + (std::rand() % (HEIGHT - 50 - MARGIN));
            const float speed = sampleMaxwellian(T, m);
            const float angle = (std::rand() % 360) / 360 * 2 * M_PI;
            vx = speed * cos(angle);
            vy = speed * sin(angle);

            validPosition = true;
            for (const auto& particle : particles) {
                Vector2D pos = particle->getPosition();
                if (Vector2D::magnitude(Vector2D(x, y) - pos) < minDistance) {
                    validPosition = false;
                    break;
                }
            }
        }

        if (ideal) {
            particles.push_back(std::make_unique<Particle>(Vector2D(x, y), Vector2D(vx, vy), m));
        } else {
            particles.push_back(std::make_unique<VWParticle>(Vector2D(x, y), Vector2D(vx, vy), m, epsilon, sigma));
        }
    }
    return particles;
}

std::vector<std::unique_ptr<Particle>> getRateParticles(Quad bounds, float dt, int rate, bool ideal, float T, float m, float epsilon, float sigma) {
    std::vector<std::unique_ptr<Particle>> particles;

    static float increment = 0;
    increment += dt * rate;
    if (increment < 1) {
        return particles;
    }





    for (int i = 0; i < round(increment); i++) {
        double x, y, vx, vy;
        x = bounds.x + (std::rand() % static_cast<int>(bounds.width));
        y = bounds.y + (std::rand() % static_cast<int>(bounds.height));

        const float speed = sampleMaxwellian(T, m);
        // Between -pi/2 and pi/2
        const float angle = (std::rand() % 180) * M_PI / 180.0 - M_PI / 2;
        vx = speed * cos(angle);
        vy = speed * sin(angle);

        if (ideal) {
            particles.push_back(std::make_unique<Particle>(Vector2D(x, y), Vector2D(vx, vy), m));
        } else {
            particles.push_back(std::make_unique<VWParticle>(Vector2D(x, y), Vector2D(vx, vy), m, epsilon, sigma));
        }
    }

    increment -= round(increment);

    return particles;
}




void createWalls(Wall walls[]) {



    // Clockwise starting from left wall
    walls[0] = Wall(
        Vector2D(MARGIN, MARGIN),
        Vector2D(MARGIN,  HEIGHT / 2 - SLIT_WIDTH / 2)
        );
    walls[1] = Wall(
        Vector2D(MARGIN, HEIGHT / 2 + SLIT_WIDTH / 2),
        Vector2D(MARGIN, HEIGHT - MARGIN)
        );
    walls[2] = Wall(
        Vector2D(MARGIN, HEIGHT - MARGIN),
        Vector2D(WIDTH / 2 - THROTTLE_LENGTH / 2, HEIGHT - MARGIN)
        );
    walls[3] = Wall(
        Vector2D(WIDTH / 2 - THROTTLE_LENGTH / 2, HEIGHT - MARGIN),
        Vector2D(WIDTH / 2 - THROTTLE_LENGTH / 2, HEIGHT - MARGIN - THROTTLE_WIDTH)
        );
    walls[4] = Wall(
        Vector2D(WIDTH / 2 - THROTTLE_LENGTH / 2, HEIGHT - MARGIN - THROTTLE_WIDTH),
        Vector2D(WIDTH / 2 + THROTTLE_LENGTH / 2, HEIGHT - MARGIN - THROTTLE_WIDTH)
        );
    walls[5] = Wall(
        Vector2D(WIDTH / 2 + THROTTLE_LENGTH / 2, HEIGHT - MARGIN - THROTTLE_WIDTH),
        Vector2D(WIDTH / 2 + THROTTLE_LENGTH / 2, HEIGHT - MARGIN)
        );
    walls[6] = Wall(
        Vector2D(WIDTH / 2 + THROTTLE_LENGTH / 2, HEIGHT - MARGIN),
        Vector2D(WIDTH - MARGIN, HEIGHT - MARGIN)
        );
    walls[7] = Wall(
        Vector2D(WIDTH - MARGIN, HEIGHT - MARGIN),
        Vector2D(WIDTH - MARGIN, HEIGHT / 2 + SLIT_WIDTH / 2)
        );
    walls[8] = Wall(
        Vector2D(WIDTH - MARGIN, HEIGHT / 2 - SLIT_WIDTH / 2),
        Vector2D(WIDTH - MARGIN, MARGIN)
        );
    walls[9] = Wall(
        Vector2D(WIDTH - MARGIN, MARGIN),
        Vector2D(WIDTH / 2 + THROTTLE_LENGTH / 2, MARGIN)
        );
    walls[10] = Wall(
        Vector2D(WIDTH / 2 + THROTTLE_LENGTH / 2, MARGIN),
        Vector2D(WIDTH / 2 + THROTTLE_LENGTH / 2, MARGIN + THROTTLE_WIDTH)
        );
    walls[11] = Wall(
        Vector2D(WIDTH / 2 + THROTTLE_LENGTH / 2, MARGIN + THROTTLE_WIDTH),
        Vector2D(WIDTH / 2 - THROTTLE_LENGTH / 2, MARGIN + THROTTLE_WIDTH)
        );
    walls[12] = Wall(
        Vector2D(WIDTH / 2 - THROTTLE_LENGTH / 2, MARGIN + THROTTLE_WIDTH),
        Vector2D(WIDTH / 2 - THROTTLE_LENGTH / 2, MARGIN)
        );
    walls[13] = Wall(
        Vector2D(WIDTH / 2 - THROTTLE_LENGTH / 2, MARGIN),
        Vector2D(MARGIN, MARGIN)
        );
}

    // Function to draw walls in the SFML window
void drawWalls(sf::RenderWindow& window, const Wall walls[], int numWalls) {
    for (int i = 0; i < numWalls; ++i) {
        sf::Vertex line[] = {
            sf::Vertex(sf::Vector2f(walls[i].getStart().x, walls[i].getStart().y)),
            sf::Vertex(sf::Vector2f(walls[i].getEnd().x, walls[i].getEnd().y))
        };
        window.draw(line, 2, sf::Lines);
    }
}


// Function to save kinetic energy data to a CSV file
void saveData(const double lT[], const double lP[], const double rT[], const double rP[], const double eT[], const int iters, const float dt) {
    std::filesystem::path currentPath = std::filesystem::current_path();
    std::ofstream file(currentPath.string() + "/ChamberTemps.csv");

    if (file.is_open()) {
        std::cout << "File opened successfully." << std::endl;
        file << "time,LT,LP,RT,RP,ET\n";
        for (int i = 0; i < iters; i++) {
            file << dt*i << "," << lT[i] << "," << lP[i] << "," << rT[i] << "," << rP[i] << "," << eT[i] <<"\n";
        }
        file.flush();
        file.close();
        std::cout << "File written and closed successfully." << std::endl;
    } else {
        std::cerr << "Unable to open file for writing: " << std::strerror(errno) << std::endl;
    }
}