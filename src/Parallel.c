#include "WCS.h"
#include "DataCube.h"
#include "Source.h"

#define BOXCAR_MIN_ITER 3
#define BOXCAR_MAX_ITER 6

CLASS DataCube
{
        char   *data;
        size_t  data_size;
        Header *header;
        int     data_type;
        int     word_size;
        size_t  dimension;
        size_t  axis_size[4];
        double  bscale;
        double  bzero;
        bool    verbosity;
};


/*PUBLIC void DataCube_boxcar_filter_1(DataCube *self,size_t radius)
{
	boxcar3d_flt_(&self->data[0],&self->axis_size[0],&self->axis_size[1],&self->axis_size[2],&radius);
}
*/
PUBLIC DataCube *fortran_DataCube_copy(const DataCube *source)
{
	DataCube *self = DataCube_new(source->verbosity);

		// Copy header
	self->header = Header_copy(source->header);

	self->data = (char *)memory(MALLOC, source->data_size, source->word_size * sizeof(char));


	// Copy remaining properties
	self->data_size    = source->data_size;
	self->data_type    = source->data_type;
	self->word_size    = source->word_size;
	self->dimension    = source->dimension;
	self->axis_size[0] = source->axis_size[0];
	self->axis_size[1] = source->axis_size[1];
	self->axis_size[2] = source->axis_size[2];
	self->axis_size[3] = source->axis_size[3];

	//Parallel data copy
	fortran_copy_data_(&source->data[0],&self->data[0],&self->axis_size[0],&self->axis_size[1],&self->axis_size[2]);
	
	return self;
}


PUBLIC void fortran_DataCube_mask_32(const DataCube *self, DataCube *maskCube, const double threshold, const int32_t value)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(maskCube);
	check_null(maskCube->data);
	ensure(self->data_type == -32 || self->data_type == -64, "Data cube must be of floating-point type.");
	ensure(maskCube->data_type == 32, "Mask cube must be of 32-bit integer type.");
	ensure(self->axis_size[0] == maskCube->axis_size[0] && self->axis_size[1] == maskCube->axis_size[1] && self->axis_size[2] == maskCube->axis_size[2], "Data cube and mask cube have different sizes.");
	ensure(threshold > 0.0, "Threshold must be positive.");
	float  thresh = threshold;

	fortran_mask32_(&self->data[0],&maskCube->data[0],&self->axis_size[0],&self->axis_size[1],&self->axis_size[2],&thresh,&value);

	return;
}
PUBLIC void fortran_DataCube_copy_blanked(DataCube *self, const DataCube *source)
{
	// Sanity checks
	check_null(self);
	check_null(self->data);
	check_null(source);
	check_null(source->data);
	ensure((self->data_type == -32 || self->data_type == -64) && (source->data_type == -32 || source->data_type == -64), "Cannot copy blanked pixels; both data cubes must be floating-point.");
	ensure(self->axis_size[0] == source->axis_size[0] && self->axis_size[1] == source->axis_size[1] && self->axis_size[2] == source->axis_size[2], "Cannot copy blanked pixels; data cubes differ in size.");
	
	// Loop over entire array and copy blanks

        //Parallel data copy
        fortran_copy_blank_(&source->data[0],&self->data[0],&self->axis_size[0],&self->axis_size[1],&self->axis_size[2]);

        return;
}


PUBLIC void fortran_DataCube_set_masked_32(DataCube *self, const DataCube *maskCube, const double value)
{
	check_null(self);
	check_null(self->data);
	check_null(maskCube);
	check_null(maskCube->data);
	ensure(self->data_type == -32 || self->data_type == -64, "Data cube must be of floating-point type.");
	ensure(maskCube->data_type == 8 || maskCube->data_type == 16 || maskCube->data_type == 32 || maskCube->data_type == 64, "Mask cube must be of integer type.");
	ensure(self->axis_size[0] == maskCube->axis_size[0] && self->axis_size[1] == maskCube->axis_size[1] && self->axis_size[2] == maskCube->axis_size[2], "Data cube and mask cube have different sizes.");
	float val_flt=value;
	fortran_sign_copy_(&maskCube->data[0],&self->data[0],&self->axis_size[0],&self->axis_size[1],&self->axis_size[2],&val_flt);

}


PUBLIC void fortran_DataCube_gaussian_filter(DataCube *self, const double sigma)
{
	long int i=0;
	int radius,niter;
	calc_filter(sigma, &radius, &niter);
		
	for(i=0;i<niter;i++)
	{
        	boxcar_x_(&self->data[0],&self->axis_size[0],&self->axis_size[1],&self->axis_size[2],&radius);
	}

	for(i=0;i<niter;i++)
	{
                boxcar_y_(&self->data[0],&self->axis_size[0],&self->axis_size[1],&self->axis_size[2],&radius);
        }
}

PUBLIC void fortran_DataCube_boxcar_filter(DataCube *self, size_t radius)
{
        long int i=0;
        boxcar_z_(&self->data[0],&self->axis_size[0],&self->axis_size[1],&self->axis_size[2],&radius);
}



/* function from SoFia2*/
void calc_filter(const double sigma, size_t *filter_radius, size_t *n_iter)
{
        *n_iter = 0;
        *filter_radius = 0;
        double tmp = -1.0;
        size_t i;

        for(i = BOXCAR_MIN_ITER; i <= BOXCAR_MAX_ITER; ++i)
        {
                const double radius = sqrt((3.0 * sigma * sigma / i) + 0.25) - 0.5;
                const double diff = fabs(radius - floor(radius + 0.5));

                if(tmp < 0.0 || diff < tmp)
                {
                        tmp = diff;
                        *n_iter = i;
                        *filter_radius = (size_t)(radius + 0.5);
                }
        }

        // Print some information
        /*const double sigma_approx = sqrt((double)(*n_iter) * ((2.0 * (double)(*filter_radius) + 1.0) * (2.0 * (double)(*filter_radius) + 1.0) - 1.0) / 12.0);
        message("Requested filter size:    sigma = %.2f\n", sigma);
        message("Approximated filter size: sigma = %.2f\n", sigma_approx);
        message("  using N = %zu and R = %zu\n", *n_iter, *filter_radius);*/

        return;
}

