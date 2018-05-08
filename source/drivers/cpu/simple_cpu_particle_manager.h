#ifndef __SIMPLE_CPU_PARTICLE_MANAGER_H_
#define __SIMPLE_CPU_PARTICLE_MANAGER_H_

#include "cpu_particle_manager.h"

namespace nbl { namespace drivers {

template<typename material_manager_t>
using simple_cpu_particle_manager = cpu_particle_manager<material_manager_t, void>;

}} // namespace nbl::drivers

#endif // __SIMPLE_CPU_PARTICLE_MANAGER_H_
