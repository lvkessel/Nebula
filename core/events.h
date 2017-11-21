#ifndef __EVENTS_H_
#define __EVENTS_H_

/*
 * This file defines possible events.
 * 
 * scatter_event represents a scattering event in the bulk of a material.
 * intersect_event represents an intersection with geometry.
 */

#include "triangle.h"

struct scatter_event
{
	uint8_t type;
	real distance;
};

struct intersect_event
{
	real isect_distance;
	triangle* isect_triangle;
};

#endif // __EVENTS_H_