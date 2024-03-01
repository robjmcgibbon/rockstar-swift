/*
 * SWIFT I/O for Rockstar
 * Matthieu Schaller (mschaller@lorentz.leidenuniv.nl)
 */
#ifdef ENABLE_HDF5

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <hdf5.h> /* HDF5 required */
#include "io_hdf5.h"
#include "io_swift.h"
#include "io_util.h"
#include "../universal_constants.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"

#define SWIFT_NTYPES 7

void swift_calculate_internal_energy(hid_t HDF_FileID, char *filename, char *gid, char *dataid, struct particle *p, int64_t to_read, int64_t offset, int64_t stride, hid_t type) {
  int64_t width = 4;
  hid_t HDF_GroupID = check_H5Gopen(HDF_FileID, gid, filename);

  // Start reading in densities
  void *density_buffer = check_malloc_s(NULL, to_read, width*stride);
  float *fdensity_buffer = density_buffer;

  hid_t HDF_DatasetID_density = check_H5Dopen(HDF_GroupID, "Densities", gid, filename);
  hid_t HDF_DataspaceID_density = check_H5Dget_space(HDF_DatasetID_density);

  check_H5Sselect_all(HDF_DataspaceID_density);
  hssize_t npoints_density = H5Sget_select_npoints(HDF_DataspaceID_density);

  if (npoints_density != to_read*stride) {
    fprintf(stderr, "[Error] dataspace %s/%s in HDF5 file %s not expected size!\n  (Actual size = %"PRId64" elements; expected size = %"PRId64" elements\n",
	    gid, "Densities", filename, (int64_t)(npoints_density), stride*to_read);
    exit(1);
  }

  check_H5Dread(HDF_DatasetID_density, type, density_buffer, "Densities", gid, filename);

  H5Sclose(HDF_DataspaceID_density);
  H5Dclose(HDF_DatasetID_density);
  // Done reading in densities

  // Start reading in pressures
  void *pressure_buffer = check_malloc_s(NULL, to_read, width*stride);
  float *fpressure_buffer = pressure_buffer;

  hid_t HDF_DatasetID_pressure = check_H5Dopen(HDF_GroupID, "Pressures", gid, filename);
  hid_t HDF_DataspaceID_pressure = check_H5Dget_space(HDF_DatasetID_pressure);

  check_H5Sselect_all(HDF_DataspaceID_pressure);
  hssize_t npoints_pressure = H5Sget_select_npoints(HDF_DataspaceID_pressure);

  if (npoints_pressure != to_read*stride) {
    fprintf(stderr, "[Error] dataspace %s/%s in HDF5 file %s not expected size!\n  (Actual size = %"PRId64" elements; expected size = %"PRId64" elements\n",
	    gid, "Pressures", filename, (int64_t)(npoints_pressure), stride*to_read);
    exit(1);
  }

  check_H5Dread(HDF_DatasetID_pressure, type, pressure_buffer, "Pressures", gid, filename);

  H5Sclose(HDF_DataspaceID_pressure);
  H5Dclose(HDF_DatasetID_pressure);
  // Done reading in pressures

  H5Gclose(HDF_GroupID);

  for (int64_t i=0; i<to_read; i++) {
    float energy = fpressure_buffer[i*stride] / (fdensity_buffer[i*stride] * 2.0 / 3.0);
    memcpy(((char *)&(p[i]))+offset, &energy, stride*width);
  }

  free(density_buffer);
  free(pressure_buffer);
}


void swift_read_dataset(hid_t HDF_FileID, char *filename, char *gid, char *dataid, struct particle *p, int64_t to_read, int64_t offset, int64_t stride, hid_t type) {
  int64_t width = (type == H5T_NATIVE_LLONG) ? 8 : 4;
  void *buffer = check_malloc_s(buffer, to_read, width*stride);
  int64_t *ibuffer = buffer;
  float *fbuffer = buffer;

  hid_t HDF_GroupID = check_H5Gopen(HDF_FileID, gid, filename);
  hid_t HDF_DatasetID = check_H5Dopen(HDF_GroupID, dataid, gid, filename);
  hid_t HDF_DataspaceID = check_H5Dget_space(HDF_DatasetID);

  check_H5Sselect_all(HDF_DataspaceID);
  hssize_t npoints = H5Sget_select_npoints(HDF_DataspaceID);

  if (npoints != to_read*stride) {
    fprintf(stderr, "[Error] dataspace %s/%s in HDF5 file %s not expected size!\n  (Actual size = %"PRId64" elements; expected size = %"PRId64" elements\n", 
	    gid, dataid, filename, (int64_t)(npoints), stride*to_read);
    exit(1);
  }

  check_H5Dread(HDF_DatasetID, type, buffer, dataid, gid, filename);
  
  H5Sclose(HDF_DataspaceID);
  H5Dclose(HDF_DatasetID);
  H5Gclose(HDF_GroupID);

  if (width == 8)
    for (int64_t i=0; i<to_read; i++)
      p[i].id = ibuffer[i];
  else
    for (int64_t i=0; i<to_read; i++)
      memcpy(((char *)&(p[i]))+offset, fbuffer+(i*stride), stride*width);

  free(buffer);
}

float swift_readattr_float(hid_t HDF_GroupID, char *filename, char* groupname, char *objName)
{
  hid_t HDF_Type = H5T_NATIVE_FLOAT;
  hid_t HDF_AttrID = check_H5Aopen_name(HDF_GroupID, objName, groupname, filename);
  hid_t HDF_DataspaceID = check_H5Aget_space(HDF_AttrID);

  check_H5Sselect_all(HDF_DataspaceID);
  
  float data = 0.0;
  check_H5Aread( HDF_AttrID, HDF_Type, &data, objName, groupname, filename);

  H5Sclose(HDF_DataspaceID);
  H5Aclose(HDF_AttrID);
  return data;
}

void swift_readattr_array(hid_t HDF_GroupID, char *filename, char *objName, hid_t type, char* groupname, void *data)
{
  hid_t HDF_AttrID = check_H5Aopen_name(HDF_GroupID, objName, groupname, filename);
  hid_t HDF_DataspaceID = check_H5Aget_space(HDF_AttrID);
  check_H5Sselect_all(HDF_DataspaceID);

  int64_t ndims = check_H5Sget_simple_extent_ndims( HDF_DataspaceID );
  assert(ndims == 1);
  hsize_t dimsize = 0;
  check_H5Sget_simple_extent_dims(HDF_DataspaceID, &dimsize);
  assert(dimsize == SWIFT_NTYPES);
  
  check_H5Aread(HDF_AttrID, type, data, objName, groupname, filename);

  H5Sclose(HDF_DataspaceID);
  H5Aclose(HDF_AttrID);
}

void swift_rescale_particles(struct particle *p, int64_t p_start, int64_t nelems) {


  // MATTHIEU: Deal with this!
  //   double vel_rescale = sqrt(SCALE_NOW)*SWIFT_VELOCITY_CONVERSION;
  // if (LIGHTCONE) vel_rescale = 1;

  // MATTHIEU: Deal with this better: Read the units from the snaps
  /* Convert to:
   * -   Mpc / h for positions
   * -   km/s for velocities
   * -   Msun / h for masses
   */
  for (int64_t i=0; i<nelems; i++) {
    for (int64_t j=0; j<3; j++) {
      p[p_start+i].pos[j]   *= SWIFT_LENGTH_CONVERSION * h0;
      /* No velocity rescale required: SWIFT vels are already a * r_dot [km/s]*/
    }
    p[p_start+i].mass *= SWIFT_MASS_CONVERSION * h0;
    p[p_start+i].energy /= (SCALE_NOW * SCALE_NOW);
  }
}

void load_particles_swift(char *filename, struct particle **p, int64_t *num_p)
{	
  hid_t HDF_FileID = check_H5Fopen(filename, H5F_ACC_RDONLY);
  hid_t HDF_Header = check_H5Gopen(HDF_FileID, "Header", filename);
  hid_t HDF_Cosmo = check_H5Gopen(HDF_FileID, "Cosmology", filename);
  
  Ol = swift_readattr_float(HDF_Cosmo, filename, "Cosmology", "Omega_lambda");
  Om = swift_readattr_float(HDF_Cosmo, filename, "Cosmology", "Omega_m");
  h0 = swift_readattr_float(HDF_Cosmo, filename, "Cosmology", "h");

  SCALE_NOW = swift_readattr_float(HDF_Header, filename, "Header", "Scale-factor");
  BOX_SIZE = swift_readattr_float(HDF_Header, filename, "Header", "BoxSize");
  
  uint32_t npart_low[SWIFT_NTYPES], npart_high[SWIFT_NTYPES] = {0};
  int64_t npart[SWIFT_NTYPES];
  
  swift_readattr_array(HDF_Header, filename, "NumPart_ThisFile", H5T_NATIVE_UINT64, "Header", npart);
  swift_readattr_array(HDF_Header, filename, "NumPart_Total_HighWord", H5T_NATIVE_UINT32, "Header", npart_high);
  swift_readattr_array(HDF_Header, filename, "NumPart_Total", H5T_NATIVE_UINT32, "Header", npart_low);

  // Convert to Msun / h 
  BOX_SIZE *= SWIFT_LENGTH_CONVERSION * h0;  
  
  H5Gclose(HDF_Header);
  H5Gclose(HDF_Cosmo);
     
  int64_t to_read = 0;
  TOTAL_PARTICLES = 0;
  for (int64_t i=0; i<SWIFT_NTYPES; i++) {
    to_read += npart[i];
    TOTAL_PARTICLES += ( ((int64_t)npart_high[i]) << 32 ) 
      + (int64_t)npart_low[i];
  }

  int64_t total_DM_particles =  (((int64_t)npart_high[1]) << 32) + (int64_t)npart_low[1]; 
  AVG_PARTICLE_SPACING = BOX_SIZE / cbrt((double) total_DM_particles); /* Mpc / h */
  PARTICLE_MASS = Om * CRITICAL_DENSITY * pow(AVG_PARTICLE_SPACING, 3.0);

  // TODO: Only print this from one process
  printf("SWIFT: filename:           %s\n", filename);
  printf("SWIFT: box size:           %g Mpc/h\n", BOX_SIZE);
  printf("SWIFT: h0:                 %g\n", h0);
  printf("SWIFT: scale factor:       %g\n", SCALE_NOW);
  printf("SWIFT: FORCE_RES:          %g Mpc/h\n", FORCE_RES);
  printf("SWIFT: FORCE_RES_PHYS_MAX: %g Mpc/h\n", FORCE_RES_PHYS_MAX);
  printf("SWIFT: Total Part:         %" PRIu64 "\n", TOTAL_PARTICLES);
  printf("SWIFT: ThisFile Part:      %" PRIu64 "\n", to_read);
  printf("SWIFT: DMO Part Mass:      %g Msun/h (DMO-equivalent if running with hydro)\n", PARTICLE_MASS);
  printf("SWIFT: avgPartSpacing:     %g Mpc/h\n\n", AVG_PARTICLE_SPACING);
  
  if (!npart[SWIFT_DM_PARTTYPE]) {
    H5Fclose(HDF_FileID);
    printf("   SKIPPING FILE, DM PARTICLE COUNT ZERO.\n");
    return;
  }

  check_realloc_s(*p, ((*num_p)+to_read), sizeof(struct particle));
  memset((*p)+(*num_p), 0, sizeof(struct particle)*to_read);

  for (int64_t i=0; i<SWIFT_NTYPES; i++) {

    /* read IDs, pos, vel */
    char buffer[100];

    /* Re-order the particle types to match the Rockstar convention */
    int32_t type = RTYPE_DM;
    if (i==4) type = RTYPE_STAR;
    else if (i==5) type = RTYPE_BH;
    else if (i==0) type = RTYPE_GAS;
    if (!npart[i]) continue;

    /* Ignore neutrinos */
    if (i > 5) break;


    snprintf(buffer, 100, "PartType%"PRId64, i);
    swift_read_dataset(HDF_FileID, filename, buffer, "ParticleIDs", *p + (*num_p),
	   npart[i], (char *)&(p[0][0].id)-(char*)(p[0]), 1, H5T_NATIVE_LLONG);
    swift_read_dataset(HDF_FileID, filename, buffer, "Coordinates", *p + (*num_p),
	   npart[i], (char *)&(p[0][0].pos[0])-(char*)(p[0]), 3, H5T_NATIVE_FLOAT);
    swift_read_dataset(HDF_FileID, filename, buffer, "Velocities", *p + (*num_p),
	   npart[i], (char *)&(p[0][0].pos[3])-(char*)(p[0]), 3, H5T_NATIVE_FLOAT);

    if (type == RTYPE_BH) {
      swift_read_dataset(HDF_FileID, filename, buffer, "DynamicalMasses", *p + (*num_p),
			 npart[i], (char *)&(p[0][0].mass)-(char*)(p[0]), 1, H5T_NATIVE_FLOAT);
    } else {
      swift_read_dataset(HDF_FileID, filename, buffer, "Masses", *p + (*num_p),
			 npart[i], (char *)&(p[0][0].mass)-(char*)(p[0]), 1, H5T_NATIVE_FLOAT);
    }
    for (int64_t j=0; j<npart[i]; j++) p[0][(*num_p)+j].type = type;

    if (type==RTYPE_GAS) {
// TODO: What if internal energy is in file?
      swift_calculate_internal_energy(HDF_FileID, filename, buffer, "InternalEnergies", *p + (*num_p),
             npart[i], (char *)&(p[0][0].energy)-(char*)(p[0]), 1, H5T_NATIVE_FLOAT);
    }

    /* if mass table is 0 but type is (primary) dark matter, need to set dark-matter particle mass */
    /* if (i == 1) { */
      
    /*   PARTICLE_MASS = p[0][*num_p].mass * SWIFT_MASS_CONVERSION; */
    /*   AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY)); */
    /* } */


    swift_rescale_particles(*p, *num_p, npart[i]);
    printf("SWIFT: first part, type %"PRId64": (%f, %f, %f, %f, %f, %f); t=%d, u=%f, m=%e\n", i, p[0][*num_p].pos[0], p[0][*num_p].pos[1], p[0][*num_p].pos[2], p[0][*num_p].pos[3], p[0][*num_p].pos[4], p[0][*num_p].pos[5], type, p[0][*num_p].energy, p[0][*num_p].mass);
    *num_p += npart[i];
    //printf("SWIFT: Particle Mass:   %g Msun/h\n", PARTICLE_MASS);
    //printf("SWIFT: avgPartSpacing: %g Mpc/h\n\n", AVG_PARTICLE_SPACING);
  }

  H5Fclose(HDF_FileID);
}

#endif /* ENABLE_HDF5 */
