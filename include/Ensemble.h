// Ensemble.h
#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <memory>
#include <SFML/Graphics/RenderWindow.hpp>

#include "Particle.h"
#include "Quad.h"
#include "Wall.h"

class Ensemble {
protected:
    std::vector<std::unique_ptr<Particle>> particles;
    Wall* bounds;
    int num_of_bounds;

public:
    Ensemble(std::vector<std::unique_ptr<Particle>> particles, Wall* bounds, int num_of_bounds);
    virtual void iterateParticles(float dt);
    void addParticles(std::vector<std::unique_ptr<Particle>> newParticles);
    [[nodiscard]] std::vector<std::unique_ptr<Particle>> const& getParticles() const { return particles; }
    [[nodiscard]] double getTemperature() const;
    [[nodiscard]] double getTemperatureInRegion(const Quad& region) const;
    [[nodiscard]] virtual double getPressureInRegion(const Quad& region) const;
    void cullNotInRegion(const Quad& region);
    void draw(sf::RenderWindow& window) const;
};

#endif // ENSEMBLE_H