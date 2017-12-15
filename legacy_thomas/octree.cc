/**
 * @file src/cdsem/octree.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#include "octree.hh"
#include <cfloat>
#include <cmath>
#include <limits>
#include <stack>
#include <algorithm>
#include "point3.hh"
#include "tribox.hh"

namespace legacy_thomas{

octree::octree(const point3& AABB_min, const point3& AABB_max) {
    _AABB_center = (AABB_max+AABB_min)/2;
    _AABB_halfsize.x = std::fabs(AABB_max.x-AABB_min.x)/2;
    _AABB_halfsize.y = std::fabs(AABB_max.y-AABB_min.y)/2;
    _AABB_halfsize.z = std::fabs(AABB_max.z-AABB_min.z)/2;
    _parent_p = nullptr;
    for(int octant = 0; octant < 8; octant++)
        _child_p[octant] = nullptr;
}

octree::octree(const octree& node)
: octree(node._AABB_center-node._AABB_halfsize, node._AABB_center+node._AABB_halfsize) {
    for(const triangle* triangle_p : node._triangles)
        insert(*triangle_p);
}

octree::~octree() {
    for(int octant = 0; octant < 8; octant++)
        if(_child_p[octant] != nullptr)
            delete _child_p[octant];
    if(_parent_p != nullptr) {
        for(int octant = 0; octant < 8; octant++)
            if(_parent_p->_child_p[octant] == this)
                _parent_p->_child_p[octant] = nullptr;
    } else {
        for(const triangle* triangle_p : _triangles)
            delete triangle_p;
    }
}

const point3& octree::center() const {
    return _AABB_center;
}

const point3& octree::halfsize() const {
    return _AABB_halfsize;
}

bool octree::is_leaf() const {
    for(int octant = 0; octant < 8; octant++)
        if(_child_p[octant] != nullptr)
            return false;
    return true;
}

int octree::level() const {
    int _level = 0;
    const octree* node_p = this;
    while(node_p->_parent_p != nullptr) {
        node_p = node_p->_parent_p;
        _level++;
    }
    return _level;
}

int octree::depth() const {
    int _depth = 0;
    std::stack<const octree*> node_p_stack;
    node_p_stack.push(this);
    while(!node_p_stack.empty()) {
        const octree* node_p = node_p_stack.top();
        node_p_stack.pop();
        _depth = std::max(_depth, node_p->level());
        for(int octant = 0; octant < 8; octant++) {
            const octree* child_p = node_p->_child_p[octant];
            if(child_p != nullptr)
                node_p_stack.push(child_p);
        }
    }
    return _depth;
}

int octree::occupancy() const {
    int _occupancy = 0;
    std::stack<const octree*> node_p_stack;
    node_p_stack.push(this);
    while(!node_p_stack.empty()) {
        const octree* node_p = node_p_stack.top();
        node_p_stack.pop();
        if(node_p->is_leaf())
            _occupancy = std::max(_occupancy, static_cast<int>(node_p->_triangles.size()));
        for(int octant = 0; octant < 8; octant++) {
            const octree* child_p = node_p->_child_p[octant];
            if(child_p != nullptr)
                if(!child_p->_triangles.empty())
                    node_p_stack.push(child_p);
        }
    }
    return _occupancy;
}

uint64_t octree::location() const {
    std::stack<int> octant_stack;
    const octree* node_p = this;
    while(node_p->_parent_p != nullptr)
        for(int octant = 0; octant < 8; octant++)
            if(node_p->_parent_p->_child_p[octant] == node_p) {
                octant_stack.push(octant);
                node_p = node_p->_parent_p;
                break;
            }
    uint64_t _location = 1;
    while(!octant_stack.empty()) {
        _location = (_location<<3)|octant_stack.top();
        octant_stack.pop();
    }
    return _location;
}

const octree* octree::root() const {
    return root();
}

octree* octree::root() {
    octree* node_p = this;
    while(node_p->_parent_p != nullptr)
        node_p = node_p->_parent_p;
    return node_p;
}

const octree* octree::parent() const {
    return _parent_p;
}

const octree* octree::traverse(int octant) const {
    return _child_p[octant%8];
}

const octree::triangle_p_vector& octree::triangles() const {
    return _triangles;
}

const triangle* octree::insert(const triangle& _triangle) {
    if(!overlap(_triangle))
        return nullptr;
    const triangle* new_triangle_p(new triangle(_triangle));
    std::stack<std::pair<octree*,const triangle*>> node_p_stack;
    node_p_stack.push(std::make_pair(root(), new_triangle_p));
    while(!node_p_stack.empty()) {
        octree* node_p = node_p_stack.top().first;
        const triangle* triangle_p = node_p_stack.top().second;
        node_p_stack.pop();
        node_p->_triangles.push_back(triangle_p);
        if(!node_p->is_leaf()) {
            /* traverse one level down */
            for(int octant = 0; octant < 8; octant++) {
                octree* child_p = node_p->_child_p[octant];
                if(child_p->overlap(*triangle_p))
                    node_p_stack.push(std::make_pair(child_p, triangle_p));
            }
        } else if((node_p->_parent_p == nullptr) || ((node_p->_triangles.size() > _split_count) && (node_p->level() < _max_depth))) {
            /* redistribute triangles to children */
            for(int octant = 0; octant < 8; octant++) {
                point3 AABB_max = node_p->_AABB_center;
                AABB_max.x += node_p->_AABB_halfsize.x*(octant&1)/1;
                AABB_max.y += node_p->_AABB_halfsize.y*(octant&2)/2;
                AABB_max.z += node_p->_AABB_halfsize.z*(octant&4)/4;
                octree* child_p = new octree(AABB_max-node_p->_AABB_halfsize, AABB_max);
                child_p->_parent_p = node_p;
                node_p->_child_p[octant] = child_p;
                for(const triangle* triangle_p : node_p->_triangles)
                    if(child_p->overlap(*triangle_p))
                        node_p_stack.push(std::make_pair(child_p, triangle_p));
            }
        }
    }
    return new_triangle_p;
}

std::pair<const triangle*,double> octree::intersect(const point3& A, const point3& B) const {
    const triangle* triangle_p = nullptr;
    double distance = 1;
    std::stack<const octree*> node_p_stack;
    if(overlap(A, B))
        node_p_stack.push(this);
    point3 C = B;
    while(!node_p_stack.empty()) {
        const octree* node_p = node_p_stack.top();
        node_p_stack.pop();
        if(node_p->is_leaf()) {
            const double eps = 10*DBL_EPSILON;
            const point3 dir = B-A;
            for(auto cit = node_p->_triangles.cbegin(); cit != node_p->_triangles.cend(); cit++) {
                const point3 e1 = (*cit)->B-(*cit)->A;
                const point3 e2 = (*cit)->C-(*cit)->A;
                // T. Möller and B. Trumbore, Journal of Graphics Tools, 2(1):21--28, 1997.
                const point3 pvec = cross_product(dir, e2);
                const double det = dot_product(e1, pvec);
                if((det > -eps) && (det < eps))
                    continue;
                const point3 tvec = A-(*cit)->A;
                const double u = dot_product(tvec, pvec)/det;
                if((u < -eps) || (u > 1.0+eps))
                    continue;
                const point3 qvec = cross_product(tvec, e1);
                const double v = dot_product(dir, qvec)/det;
                if((v < -eps) || (u+v > 1.0+eps))
                    continue;
                const double t = dot_product(e2, qvec)/det;
                if((t > 0) && (t <= distance)) {
                    triangle_p = *cit;
                    distance = t;
                }
            }
            if(triangle_p != nullptr)
                C = A+dir*distance;
        }
        for(int octant = 0; octant < 8; octant++) {
            const octree* child_p = node_p->_child_p[octant];
            if(child_p != nullptr)
                if(!child_p->_triangles.empty())
                    if(child_p->overlap(A, C))
                        node_p_stack.push(child_p);
        }
    }
    return std::make_pair(triangle_p, distance);
}

int octree::octant(const point3& pos) const {
    int octant = 0;
    octant += (pos.x > _AABB_center.x) ? 1 : 0;
    octant += (pos.y > _AABB_center.y) ? 2 : 0;
    octant += (pos.z > _AABB_center.z) ? 4 : 0;
    return octant;
}

bool octree::overlap(const triangle& _triangle) const {
    double boxcenter[3] = {_AABB_center.x, _AABB_center.y, _AABB_center.z};
    // LvK 20170420: include edge-case triangles. TODO
    double boxhalfsize[3] = {1.01*_AABB_halfsize.x, 1.01*_AABB_halfsize.y, 1.01*_AABB_halfsize.z};
	//double boxhalfsize[3] = { _AABB_halfsize.x, _AABB_halfsize.y, _AABB_halfsize.z };
    double triverts[3][3];
    triverts[0][0] = _triangle.A.x; triverts[0][1] = _triangle.A.y; triverts[0][2] = _triangle.A.z;
    triverts[1][0] = _triangle.B.x; triverts[1][1] = _triangle.B.y; triverts[1][2] = _triangle.B.z;
    triverts[2][0] = _triangle.C.x; triverts[2][1] = _triangle.C.y; triverts[2][2] = _triangle.C.z;
    if(triBoxOverlap(boxcenter, boxhalfsize, triverts) > 0)
        return true;
    return false;
}

bool octree::overlap(const point3& A, const point3& B) const {
    // A. Williams, S. Barrus, et al., Journal of Graphics Tools, 10(1):49-54, 2005
    const point3 bmin = _AABB_center-_AABB_halfsize;
    const point3 bmax = _AABB_center+_AABB_halfsize;
    const double divx = 1.0/(B.x-A.x);
    const double divy = 1.0/(B.y-A.y);
    const double divz = 1.0/(B.z-A.z);
    double tmin, tmax;
    if(divx >= 0) {
        tmin = (bmin.x-A.x)*divx;
        tmax = (bmax.x-A.x)*divx;
    } else {
        tmin = (bmax.x-A.x)*divx;
        tmax = (bmin.x-A.x)*divx;
    }
    double tymin, tymax;
    if(divy >= 0) {
        tymin = (bmin.y-A.y)*divy;
        tymax = (bmax.y-A.y)*divy;
    } else {
        tymin = (bmax.y-A.y)*divy;
        tymax = (bmin.y-A.y)*divy;
    }
    if((tmin > tymax) || (tymin > tmax))
        return false;
    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);
    double tzmin, tzmax;
    if(divz >= 0) {
        tzmin = (bmin.z-A.z)*divz;
        tzmax = (bmax.z-A.z)*divz;
    } else {
        tzmin = (bmax.z-A.z)*divz;
        tzmax = (bmin.z-A.z)*divz;
    }
    if((tmin > tzmax) || (tzmin > tmax))
        return false;
    tmin = std::max(tmin, tzmin);
    tmax = std::min(tmax, tzmax);
    if((tmax < 0) || (tmin > 1))
        return false;
    return true;
}


}
