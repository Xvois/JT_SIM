//
// Created by Sonny Parker on 25/11/2024.
//
#include "../include/QuadTree.h"

QuadTree::QuadTree(const Quad& boundary) : boundary(boundary), divided(false) {}

void QuadTree::subdivide() {
    float x = boundary.x;
    float y = boundary.y;
    float w = boundary.width / 2;
    float h = boundary.height / 2;

    northwest = std::make_unique<QuadTree>(Quad(x - w / 2, y - h / 2, w, h));
    northeast = std::make_unique<QuadTree>(Quad(x + w / 2, y - h / 2, w, h));
    southwest = std::make_unique<QuadTree>(Quad(x - w / 2, y + h / 2, w, h));
    southeast = std::make_unique<QuadTree>(Quad(x + w / 2, y + h / 2, w, h));
    divided = true;
}

void QuadTree::clear() {
    particles.clear();
    divided = false;
    northwest.reset();
    northeast.reset();
    southwest.reset();
    southeast.reset();
}

bool QuadTree::insert(Particle* particle) {
    if (!boundary.contains(particle->getPosition().x, particle->getPosition().y)) {
        return false;
    }

    if (particles.size() < capacity) {
        particles.push_back(particle);
        return true;
    } else {
        if (!divided) {
            subdivide();
        }

        if (northwest->insert(particle)) return true;
        if (northeast->insert(particle)) return true;
        if (southwest->insert(particle)) return true;
        if (southeast->insert(particle)) return true;
    }
    return false;
}

void QuadTree::query(const Quad& range, std::vector<Particle*>& found) const {
    if (!boundary.intersects(range)) {
        return;
    }

    for (Particle* p : particles) {
        if (range.contains(p->getPosition().x, p->getPosition().y)) {
            found.push_back(p);
        }
    }

    if (divided) {
        northwest->query(range, found);
        northeast->query(range, found);
        southwest->query(range, found);
        southeast->query(range, found);
    }
}
