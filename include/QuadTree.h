//
// Created by Sonny Parker on 25/11/2024.
//

#ifndef QUADTREE_H
#define QUADTREE_H

#include <vector>

#include "Particle.h"
#include "Quad.h"

class QuadTree {
private:
    static const int capacity = 4;
    Quad boundary;
    std::vector<Particle*> particles;
    bool divided;
    std::unique_ptr<QuadTree> northwest;
    std::unique_ptr<QuadTree> northeast;
    std::unique_ptr<QuadTree> southwest;
    std::unique_ptr<QuadTree> southeast;

    void subdivide();

public:
    QuadTree(const Quad& boundary);
    void clear();
    bool insert(Particle* particle);
    void query(const Quad& range, std::vector<Particle*>& found) const;
};


#endif //QUADTREE_H
