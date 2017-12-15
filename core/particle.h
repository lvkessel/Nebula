#ifndef __PARTICLE_H_
#define __PARTICLE_H_

struct particle
{
	vec3 pos;        // Current position
	vec3 dir;        // Current direction (unnormalised)
	real kin_energy; // Current kinetic energy
};

#endif // __PARTICLE_H_
