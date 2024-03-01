#ifndef _IO_SWIFT_H_
#define _IO_SWIFT_H_
#ifdef ENABLE_HDF5
#include <stdint.h>
#include "../particle.h"
#include <hdf5.h>

void load_particles_swift(char *filename, struct particle **p, int64_t *num_p);
void swift_read_dataset(hid_t HDF_FileID, char *filename, char *gid, char *dataid, struct particle *p, int64_t to_read, int64_t offset, int64_t stride, hid_t type);

#endif /* ENABLE_HDF5 */
#endif /* _IO_SWIFT_H_ */
