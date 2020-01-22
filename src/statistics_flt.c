/// ____________________________________________________________________ ///
///                                                                      ///
/// SoFiA 2.0.1 (statistics_flt.c) - Source Finding Application          ///
/// Copyright (C) 2019 Tobias Westmeier                                  ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// Address:  Tobias Westmeier                                           ///
///           ICRAR M468                                                 ///
///           The University of Western Australia                        ///
///           35 Stirling Highway                                        ///
///           Crawley WA 6009                                            ///
///           Australia                                                  ///
///                                                                      ///
/// E-mail:   tobias.westmeier [at] uwa.edu.au                           ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// This program is free software: you can redistribute it and/or modify ///
/// it under the terms of the GNU General Public License as published by ///
/// the Free Software Foundation, either version 3 of the License, or    ///
/// (at your option) any later version.                                  ///
///                                                                      ///
/// This program is distributed in the hope that it will be useful,      ///
/// but WITHOUT ANY WARRANTY; without even the implied warranty of       ///
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ///
/// GNU General Public License for more details.                         ///
///                                                                      ///
/// You should have received a copy of the GNU General Public License    ///
/// along with this program. If not, see http://www.gnu.org/licenses/.   ///
/// ____________________________________________________________________ ///
///                                                                      ///

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "statistics_flt.h"



// -------------------------------- //
// Check if data array contains NaN //
// -------------------------------- //

int contains_nan_flt(const float *data, const size_t size)
{
	const float *ptr = data + size;
	while(ptr --> data) if(IS_NAN(*ptr)) return 1;
	return 0;
}



// --------------------------------------------------------- //
// Maximum and minimum                                       //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1)      data - Pointer to the input data array         //
//   (2)      size - Size of the input data array            //
//   (3) value_max - Pointer to variable for holding maximum //
//   (4) value_min - Pointer to variable for holding minimum //
//                                                           //
// Returns:                                                  //
//                                                           //
//   No return value.                                        //
//                                                           //
// Description:                                              //
//                                                           //
//   Determines the maximum and minimum value in the input   //
//   data array simultaneously. This is faster that calling  //
//   the max() and min() functions separately in situations  //
//   where both values are required. If the array contains   //
//   only NaN, value_min and value_max are both set to NaN.  //
// --------------------------------------------------------- //

void max_min_flt(const float *data, const size_t size, float *value_max, float *value_min)
{
	const float *tmp1;
	const float *tmp2;
	const float *data2 = data + 1;
	const float *ptr = data + size - 1;
	const float *result_min;
	const float *result_max;
	
	// Find last non-NaN element
	while(IS_NAN(*ptr) && data <-- ptr);
	result_max = ptr;
	
	if(IS_NAN(*result_max))
	{
		// Data array only contains NaN
		*value_max = NAN;
		*value_min = NAN;
		return;
	}
	
	if(ptr == data)
	{
		// Data array only contains one non-NaN value
		*value_max = *result_max;
		*value_min = *result_max;
		return;
	}
	
	// Find second-to-last non-NaN element
	if(--ptr > data) while(IS_NAN(*ptr) && data <-- ptr);
	result_min = ptr;
	
	if(IS_NAN(*result_min))
	{
		// Data array only contains one non-NaN value
		*value_max = *result_max;
		*value_min = *result_max;
		return;
	}
	
	// Swap min and max if necessary
	if(*result_min > *result_max)
	{
		tmp1 = result_min;
		result_min = result_max;
		result_max = tmp1;
	}
	
	// Iterate over remaining array to update min and max
	while(ptr > data2)
	{
		tmp1 = --ptr;
		tmp2 = --ptr;
		
		// Bizarrely, the following faster by more than a factor of 2
		// than first checking if tmp1 > tmp2, presumably due to branch
		// predictor optimisation.
		if(*tmp1 > *result_max || *tmp2 > *result_max)
		{
			if(*tmp2 > *tmp1) result_max = tmp2;
			else result_max = tmp1;
		}
		if(*tmp1 < *result_min || *tmp2 < *result_min)
		{
			if(*tmp2 < *tmp1) result_min = tmp2;
			else result_min = tmp1;
		}
	}
	
	// Ensure that we didn't miss the first array element
	if(ptr > data)
	{
		tmp1 = --ptr;
		if(*tmp1 > *result_max) result_max = tmp1;
		if(*tmp1 < *result_min) result_min = tmp1;
	}
	
	*value_max = *result_max;
	*value_min = *result_min;
	
	return;
}



// --------------------------------------------------------- //
// Maximum                                                   //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the input data array              //
//   (2) size - Size of the input data array                 //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Maximum value of the input array. NaN will be returned  //
//   for NaN-only input arrays.                              //
//                                                           //
// Description:                                              //
//                                                           //
//   Determines the maximum value in the input data array.   //
//   If the array contains only NaN, then NaN is returned.   //
// --------------------------------------------------------- //

float max_flt(const float *data, const size_t size)
{
	const float *result = data + size - 1;
	while(IS_NAN(*result) && data <-- result);
	const float *ptr = result;
	
	while(ptr --> data)
	{
		// Note that this weird construct is a lot faster than without
		// the useless if clause, possibly due to the branch predictor.
		if(*ptr < *result);
		else if(*ptr > *result) result = ptr;
	}
	
	return *result;
}



// --------------------------------------------------------- //
// Minimum                                                   //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the input data array              //
//   (2) size - Size of the input data array                 //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Minimum value of the input array. NaN will be returned  //
//   for NaN-only input arrays.                              //
//                                                           //
// Description:                                              //
//                                                           //
//   Determines the minimum value in the input data array.   //
//   If the array contains only NaN, then NaN is returned.   //
// --------------------------------------------------------- //

float min_flt(const float *data, const size_t size)
{
	const float *result = data + size - 1;
	while(IS_NAN(*result) && data <-- result);
	const float *ptr = result;
	
	while(ptr --> data)
	{
		// Note that this weird construct is a lot faster than without
		// the useless if clause, possibly due to the branch predictor.
		if(*ptr > *result);
		else if(*ptr < *result) result = ptr;
	}
	
	return *result;
}



// --------------------------------------------------------- //
// Sum and mean                                              //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the input data array              //
//   (2) size - Size of the input data array                 //
//   (3) mean - If true, return the mean; otherwise the sum  //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Sum or mean of the data array. NaN will be returned for //
//   NaN-only input arrays.                                  //
//                                                           //
// Description:                                              //
//                                                           //
//   Depending on the value of mean, either the sum or the   //
//   mean of the non-NaN elements of the input array is      //
//   returned. Two convenient wrapper functions, sum() and   //
//   mean(), have been provided to explicitly return the sum //
//   and the mean of the array, respectively.                //
// --------------------------------------------------------- //

double summation_flt(const float *data, const size_t size, const bool mean)
{
	double result = 0.0;
	size_t counter = 0;
	const float *ptr = data + size;
	
	while(ptr --> data)
	{
		if(IS_NOT_NAN(*ptr))
		{
			result += *ptr;
			++counter;
		}
	}
	
	if(counter) return mean ? result / counter : result;
	return NAN;
}
double sum_flt(const float *data, const size_t size) { return summation_flt(data, size, false); }
double mean_flt(const float *data, const size_t size) { return summation_flt(data, size, true); }



// --------------------------------------------------------- //
// N-th moment about value                                   //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1)  data - Pointer to the input data array             //
//   (2)  size - Size of the input data array                //
//   (3) order - Order of the moment to be calculated        //
//   (4) value - Value about which to calculate the moment   //
//                                                           //
// Returns:                                                  //
//                                                           //
//   N-th moment of the data array. NaN will be returned for //
//   NaN-only input arrays.                                  //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the N-th moment of the input data array. The //
//   N-th moment is defined as                               //
//                                                           //
//     sum((x - value)^N) / n                                //
//                                                           //
//   where the summation is over all n data values, x, of    //
//   the input array that are not NaN. The most sensible use //
//   would be to calculate the moment about the mean.        //
// --------------------------------------------------------- //

double moment_flt(const float *data, const size_t size, unsigned int order, const double value)
{
	if(order == 0) return 1.0;
	
	const float *ptr = data + size;
	double result = 0.0;
	size_t counter = 0;
	double tmp;
	unsigned int i;
	
	while(ptr --> data)
	{
		if(IS_NOT_NAN(*ptr))
		{
			tmp = *ptr - value;
			for(i = order; --i;) tmp *= (*ptr - value);
			result += tmp;
			++counter;
		}
	}
	
	if(counter) return result / counter;
	else return NAN;
}



// --------------------------------------------------------- //
// Moments 2, 3 and 4 about value                            //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1)  data - Pointer to the input data array             //
//   (2)  size - Size of the input data array                //
//   (3) value - Value about which to calculate the moments  //
//   (4)    m2 - Pointer to a variable for storing moment 2  //
//   (5)    m3 - Pointer to a variable for storing moment 3  //
//   (6)    m4 - Pointer to a variable for storing moment 4  //
//                                                           //
// Returns:                                                  //
//                                                           //
//   No return value.                                        //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the 2nd, 3rd and 4th moment of a data array  //
//   simultaneously. This will enable faster calculation of  //
//   skewness and kurtosis, both of which depend on multiple //
//   moments of the same data array. Moments are defined in  //
//   the same way as in function moment().                   //
// --------------------------------------------------------- //

void moments_flt(const float *data, const size_t size, const double value, double *m2, double *m3, double *m4)
{
	const float *ptr = data + size;
	double result2 = 0.0;
	double result3 = 0.0;
	double result4 = 0.0;
	size_t counter = 0;
	double tmp, tmp2;
	
	while(ptr --> data)
	{
		if(IS_NOT_NAN(*ptr))
		{
			tmp = *ptr - value;
			tmp2 = tmp * tmp;
			result2 += tmp2;
			result3 += tmp2 * tmp;
			result4 += tmp2 * tmp2;
			++counter;
		}
	}
	
	if(counter)
	{
		*m2 = result2 / counter;
		*m3 = result3 / counter;
		*m4 = result4 / counter;
	}
	else
	{
		*m2 = NAN;
		*m3 = NAN;
		*m4 = NAN;
	}
	
	return;
}



// --------------------------------------------------------- //
// Standard deviation about value                            //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1)    data - Pointer to the input data array           //
//   (2)    size - Size of the input data array              //
//   (3)   value - Value about which to calculate the        //
//                 standard deviation                        //
//   (4) cadence - Can be set to > 1 to speed up algorithm   //
//   (4)   range - Flux range to be used. Can be Negative    //
//                 (-1), Full (0) or Positive (1).           //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Standard deviation about the specified value of the     //
//   elements in the input data array. NaN will be returned  //
//   for NaN-only input arrays.                              //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the standard deviation about a user-speci-   //
//   fied value of all elements in the input data array that //
//   are not NaN. The standard deviation is defined as       //
//                                                           //
//     sqrt[sum(x - value)^2 / n]                            //
//                                                           //
//   where the summation is over all n non-NaN elements, x,  //
//   of the input array.                                     //
// --------------------------------------------------------- //

double std_dev_val_flt(const float *data, const size_t size, const double value, const size_t cadence, const int range)
{
	const float *ptr = data + size;
	const float *ptr2 = data + cadence - 1;
	double result = 0.0;
	size_t counter = 0;
	double tmp;
	
	while(ptr > ptr2)
	{
		ptr -= cadence;
		
		if((range == 0 && IS_NOT_NAN(*ptr)) || (range < 0 && *ptr < 0.0) || (range > 0 && *ptr > 0.0))
		{
			tmp = (*ptr - value);
			result += tmp * tmp;
			++counter;
		}
	}
	
	return counter ? sqrt(result / counter) : NAN;
}



// --------------------------------------------------------- //
// Standard deviation about mean                             //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the input data array              //
//   (2) size - Size of the input data array                 //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Standard deviation about the mean of the elements in    //
//   the input data array. NaN will be returned for NaN-only //
//   input arrays.                                           //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the standard deviation about the mean of all //
//   elements in the input data array that are not NaN. The  //
//   standard deviation is defined as                        //
//                                                           //
//     sqrt{sum[x - mean(x)]^2 / n}                          //
//                                                           //
//   where the summation is over all n non-NaN elements, x,  //
//   of the input array.                                     //
// --------------------------------------------------------- //

double std_dev_flt(const float *data, const size_t size)
{
	return std_dev_val_flt(data, size, mean_flt(data, size), 1, 0);
}



// --------------------------------------------------------- //
// N-th smallest element in array                            //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the data array to be sorted       //
//   (2) size - Size of the input array                      //
//   (3)    n - n-th smallest value will be returned         //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Value of the n-th smallest array element.               //
//                                                           //
// Description:                                              //
//                                                           //
//   Partially sorts the input data array until the n-th-    //
//   smallest element is at the n-th position within the     //
//   array. After sorting, all elements below position n     //
//   will be smaller than or equal to the n-th-smallest      //
//   element, while  all elements above position n will be   //
//   greater than or equal to the n-th-smallest element. The //
//   value of the n-th-smallest element is returned.         //
//   Note that this function is not NaN-safe and will modify //
//   the original data array.                                //
// --------------------------------------------------------- //

float nth_element_flt(float *data, const size_t size, const size_t n)
{
	float *l = data;
	float *m = data + size - 1;
	float *ptr = data + n;
	
	while(l < m)
	{
		float value = *ptr;
		float *i = l;
		float *j = m;
		
		do
		{
			while(*i < value) ++i;
			while(value < *j) --j;
			
			if(i <= j)
			{
				float tmp = *i;
				*i = *j;
				*j = tmp;
				++i;
				--j;
			}
		} while(i <= j);
		
		if(j < ptr) l = i;
		if(ptr < i) m = j;
	}
	
	return *ptr;
}



// --------------------------------------------------------- //
// Median                                                    //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the data array to be sorted       //
//   (2) size - Size of the input array                      //
//   (3) fast - If true, will return approximate median for  //
//              even-sized arrays. If false, the exact       //
//              median will be returned, which takes longer. //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Median of the data array values.                        //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the exact median of the input data array.    //
//   NOTE that this function is not NaN-safe and will modify //
//   the original data array.                                //
// --------------------------------------------------------- //

float median_flt(float *data, const size_t size, const bool fast)
{
	const float value = nth_element_flt(data, size, size / 2);
	return IS_ODD(size) || fast ? value : (value + max_flt(data, size / 2)) / 2.0;
}

// Same, but does not alter the data array.

float median_safe_flt(float *data, const size_t size, const bool fast)
{
	float *data_copy = (float *)memory(MALLOC, size, sizeof(float));
	memcpy(data_copy, data, size * sizeof(float));
	const float result = median_flt(data_copy, size, fast);
	free(data_copy);
	return result;
}



// --------------------------------------------------------- //
// Median absolute deviation from value                      //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1)  data   - Pointer to the data array to be sorted    //
//   (2)  size   - Size of the input array                   //
//   (3) value   - Value about which to calculate the MAD    //
//   (4) cadence - Can be set to > 1 to speed up algorithm   //
//   (5)   range - Flux range to be used. Can be Negative    //
//                 (-1), Full (0) or Positive (1).           //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Median absolute deviation of the array values from the  //
//   specified data value.                                   //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the median absolute deviation (MAD) of the   //
//   data array values from a user-specified value. The MAD  //
//   is defined as                                           //
//                                                           //
//     median(|x - value|)                                   //
//                                                           //
//   where x denotes the data values from the input array.   //
//   NOTE that this function IS NaN-safe and will NOT modify //
//   the original data array.                                //
// --------------------------------------------------------- //

float mad_val_flt(const float *data, const size_t size, const float value, const size_t cadence, const int range)
{
	// Create copy of data array with specified range and cadence
	const size_t data_copy_size = (range == 0) ? (size / cadence) : (size / (2 * cadence));
	float *data_copy = (float *)memory(MALLOC, data_copy_size, sizeof(float));
	
	// Some settings
	const float *ptr = data + size;
	float *ptr_copy = data_copy;
	size_t counter = 0;
	
	// Copy |*ptr - value| into array copy
	while((ptr -= cadence) > data && counter < data_copy_size)
	{
		if((range < 0 && *ptr < 0.0) || (range == 0 && IS_NOT_NAN(*ptr)) || (range > 0 && *ptr > 0.0))
		{
			ptr_copy[counter] = fabs(*ptr - value);
			++counter;
		}
	}
	
	// Determine median
	const float result = median_flt(data_copy, counter, false);
	
	// Clean up
	free(data_copy);
	
	return result;
}



// --------------------------------------------------------- //
// Median absolute deviation                                 //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the data array to be sorted       //
//   (2) size - Size of the input array                      //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Median absolute deviation of the array values.          //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the median absolute deviation (MAD) of the   //
//   data array values. The MAD is defined as                //
//                                                           //
//     median(|x - median(x)|)                               //
//                                                           //
//   where x denotes the data values from the input array.   //
//   In the case of normally distributed random data values  //
//   the standard deviation of the data values about the     //
//   mean can be deduced by multiplying the MAD with the     //
//   constant MAD_TO_STD.                                    //
//   Note that this function is NaN-safe and will not modify //
//   the original data array.                                //
// --------------------------------------------------------- //

float mad_flt(float *data, const size_t size)
{
	return mad_val_flt(data, size, median_flt(data, size, false), 1, 0);
}



// --------------------------------------------------------- //
// Robust noise measurement for contiguous data array        //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//    (1) data - Pointer to the data array.                  //
//    (2) size - Size of the array.                          //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Robust noise measurement of the data.                   //
//                                                           //
// Description:                                              //
//                                                           //
//   Uses a robust way of measuring the noise within the     //
//   specified data array based on the median absolute de-   //
//   viation (MAD). The following assumptions and conditions //
//   apply:                                                  //
//                                                           //
//    - The noise is Gaussian and centred on zero.           //
//    - Only negative data points will be used in the noise  //
//      measurement.                                         //
//    - The algorithm is NAN-safe and discards all NANs.     //
//                                                           //
//   The result of the median of the absolute values of the  //
//   negative elements in the array will be returned. If no  //
//   valid data are found, NaN will instead be returned.     //
// --------------------------------------------------------- //

float robust_noise_flt(const float *data, const size_t size)
{
	// Allocate memory for 1D data copy
	float *data_copy = (float *)memory(MALLOC, size, sizeof(float));
	float *ptr_copy = data_copy;
	
	// Copy values of negative elements
	for(const float *ptr = data + size; ptr --> data;) if(*ptr < 0.0) *ptr_copy++ = *ptr;
	
	// Calculate median
	const size_t size_copy = ptr_copy - data_copy;
	const float result = size_copy ? -MAD_TO_STD * nth_element_flt(data_copy, size_copy, size_copy / 2) : NAN;
	// NOTE: Multiplication by -1 because values are all negative.
	
	// Clean up
	free(data_copy);
	return result;
}



// --------------------------------------------------------- //
// Robust noise measurement in region of 3D array            //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//    (1) data - Pointer to the data array.                  //
//    (2) nx   - Size of the data array in x.                //
//    (3) ny   - Size of the data array in y.                //
//    (4) nz   - Size of the data array in z.                //
//    (5) x1   - Lower boundary of region in x.              //
//    (6) x2   - Upper boundary of region in x.              //
//    (7) y1   - Lower boundary of region in y.              //
//    (8) y2   - Upper boundary of region in y.              //
//    (9) z1   - Lower boundary of region in z.              //
//   (10) z2   - Upper boundary of region in z.              //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Robust noise measurement within the specified region.   //
//                                                           //
// Description:                                              //
//                                                           //
//   Uses a robust way of measuring the noise within the     //
//   specified region of a 3D data array based on the median //
//   absolute deviation (MAD). The following assumptions     //
//   and conditions apply:                                   //
//                                                           //
//    - The array is contiguous in x and least contiguous    //
//      along the z axis.                                    //
//    - The noise is Gaussian and centred on zero.           //
//    - Only negative data points will be used in the noise  //
//      measurement.                                         //
//    - The algorithm is NAN-safe and discards all NANs.     //
//                                                           //
//   The result of the median of the absolute values of the  //
//   negative elements in the region will be returned. If no //
//   valid data are found, NaN will instead be returned.     //
// --------------------------------------------------------- //

float robust_noise_in_region_flt(const float *data, const size_t nx, const size_t ny, const size_t nz, const size_t x1, const size_t x2, const size_t y1, const size_t y2, const size_t z1, const size_t z2)
{
	// Allocate memory for 1D data copy
	float *data_copy = (float *)memory(MALLOC, (x2 - x1 + 1) * (y2 - y1 + 1) * (z2 - z1 + 1), sizeof(float));
	float *ptr = data_copy;
	float value;
	
	// Copy absolute values of negative elements
	for(size_t z = z1; z <= z2; ++z)
	{
		for(size_t y = y1; y <= y2; ++y)
		{
			for(size_t x = x1; x <= x2; ++x)
			{
				value = data[x + nx * (y + ny * z)];
				if(value < 0.0) *ptr++ = -value;
			}
		}
	}
	
	// Calculate median
	const size_t size = ptr - data_copy;
	const float result = size ? MAD_TO_STD * nth_element_flt(data_copy, size, size / 2) : NAN;
	
	// Clean up
	free(data_copy);
	return result;
}



// --------------------------------------------------------- //
// Create histogram from data array                          //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data     - Pointer to the data array to be sorted.  //
//   (2) size     - Size of the input array.                 //
//   (3) n_bins   - Number of bins of the histogram.         //
//   (4) data_min - Lower flux limit of the histogram.       //
//   (5) data_max - Upper flux limit of the histogram.       //
//   (6) cadence - Cadence for generation of histogram. A    //
//                 cadence of N means that only every N-th   //
//                 data point will be used.                  //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Pointer to the generated histogram.                     //
//                                                           //
// Description:                                              //
//                                                           //
//   Generates a histogram of the data values in 'data' with //
//   the parameters as specified by the user. A pointer to   //
//   the generated histogram will be returned. NOTE that it  //
//   is the responsibility of the user to de-allocate the    //
//   memory occupied by the histogram once it is no longer   //
//   required (using free()) to ensure that there are no me- //
//   mory leaks. This function is NaN-safe, and NaN values   //
//   will simply be ignored (due to >= and <= comparison).   //
// --------------------------------------------------------- //

size_t *create_histogram_flt(const float *data, const size_t size, const size_t n_bins, const float data_min, const float data_max, const size_t cadence)
{
	// Allocate memory
	size_t *histogram = (size_t *)memory(CALLOC, n_bins, sizeof(size_t));
	
	// Basic setup
	const float slope = (float)(n_bins - 1) / (data_max - data_min);
	const float *ptr  = data + size;
	const float tmp   = 0.5 - slope * data_min;  // The 0.5 is needed for correct rounding later on.
	
	// Generate histogram
	while((ptr -= cadence) > data)
	{
		if(*ptr >= data_min && *ptr <= data_max)
		{
			size_t bin = (size_t)(slope * *ptr + tmp);
			++histogram[bin];
		}
	}
	
	return histogram;
}



// --------------------------------------------------------- //
// Gaussian fit to histogram                                 //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data    - Pointer to the data array to be sorted.   //
//   (2) size    - Size of the input array.                  //
//   (3) cadence - Cadence for generation of histogram. A    //
//                 cadence of N means that only every N-th   //
//                 data point will be used.                  //
//   (4) range   - Flux range to use in histogram. Can be    //
//                 -1, 0 or +1 to use only negative pixels,  //
//                 all pixels, or only positive pixels, re-  //
//                 spectively.                               //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Standard deviation from Gaussian fit.                   //
//                                                           //
// Description:                                              //
//                                                           //
//   Generates a histogram from the data values in the input //
//   array and fits a Gaussian function to that histogram to //
//   determine the standard deviation. This is useful for    //
//   measuring the noise level in the data cube.             //
//                                                           //
//   The user can choose which pixels from the data cube     //
//   contribute to the histogram (positive, negative, or all //
//   pixels and every N-th pixel as controlled by cadence).  //
//   For the purpose of fitting, the histogram range will be //
//   set such that the second moment fills an optimal frac-  //
//   tion of the histogram width. A linear regression to the //
//   logarithmic histogram, ln(h(x)) = a x^2 + b, will then  //
//   be performed to derive the standard deviation from the  //
//   slope of the fit. The advantage of this technique is    //
//   that the regression can be solved analytically, making  //
//   the algorithm extremely fast. The disadvantage is that  //
//   the weighting of histogram bins will be higher near the //
//   centre of the flux distribution and lower in the wings, //
//   although this can actually be an advantage, as the      //
//   wings are often affected by non-Gaussian effects such   //
//   as flux from actual objects.                            //
// --------------------------------------------------------- //

float gaufit_flt(const float *data, const size_t size, const size_t cadence, const int range)
{
	// Determine maximum and minimum
	float data_max = 0.0;
	float data_min = 0.0;
	max_min_flt(data, size, &data_max, &data_min);
	
	if(data_min >= 0.0 || data_max <= 0.0)
	{
		warning("Maximum is not greater than minimum.");
		return NAN;
	}
	
	// Determine data range
	if(range < 0)
	{
		if(data_min >= 0.0)
		{
			warning("Minimum is not less than zero.");
			return NAN;
		}
		data_max = 0.0;
	}
	else if(range > 0)
	{
		if(data_max <= 0.0)
		{
			warning("Maximum is not greater than zero.");
			return NAN;
		}
		data_min = 0.0;
	}
	else
	{
		const float limit = fabs(data_min) < fabs(data_max) ? fabs(data_min) : fabs(data_max);
		data_min = -limit;
		data_max = limit;
	}
	
	// Create initial histogram
	const size_t n_bins = 101;
	size_t origin = n_bins / 2;
	if(range < 0) origin = n_bins - 1;
	else if(range > 0) origin = 0;
	const float inv_optimal_mom2 = 5.0 / n_bins;  // Require standard deviation to cover 1/5th of the histogram for optimal results
	size_t *histogram = create_histogram_flt(data, size, n_bins, data_min, data_max, cadence);
	
	// Calculate second moment
	float mom0 = 0.0;
	float mom1 = 0.0;
	float mom2 = 0.0;
	
	for(size_t i = n_bins; i--;)
	{
		mom0 += histogram[i];
		mom1 += histogram[i] * (float)i;
	}
	mom1 /= mom0;
	
	for(size_t i = n_bins; i--;) mom2 += histogram[i] * (mom1 - i) * (mom1 - i);
	mom2 = sqrt(mom2 / mom0);
	
	// Ensure that 2nd moment is equal to optimal_mom2 for optimal fitting
	if(range < 0) data_min *= mom2 * inv_optimal_mom2;
	else if(range > 0) data_max *= mom2 * inv_optimal_mom2;
	else
	{
		data_min *= mom2 * inv_optimal_mom2;
		data_max *= mom2 * inv_optimal_mom2;
	}
	
	// Regenerate histogram with new scaling
	free(histogram);
	histogram = create_histogram_flt(data, size, n_bins, data_min, data_max, cadence);
	
	// Fit Gaussian function using linear regression
	// (excluding first and last point in case of edge effects)
	float mean_x = 0.0;
	float mean_y = 0.0;
	size_t counter = 0;
	
	for(size_t i = n_bins - 1; i --> 1;)
	{
		if(histogram[i])
		{
			long int ii = i - origin;
			mean_x += (float)(ii * ii);
			mean_y += log((float)(histogram[i]));
			++counter;
		}
	}
	
	mean_x /= counter;
	mean_y /= counter;
	
	float upper_sum = 0.0;
	float lower_sum = 0.0;
	
	for(size_t i = n_bins - 1; i --> 1;)
	{
		if(histogram[i])
		{
			long int ii = i - origin;
			const float x = (float)(ii * ii);
			const float y = log((float)(histogram[i]));
			upper_sum += (x - mean_x) * (y - mean_y);
			lower_sum += (x - mean_x) * (x - mean_x);
		}
	}
	
	// Determine standard deviation from slope
	const float sigma = sqrt(-0.5 * lower_sum / upper_sum) * (data_max - data_min) / (n_bins - 1);
	
	// Clean up
	free(histogram);
	
	return sigma;
}



// --------------------------------------------------------- //
// Skewness                                                  //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the data array to be sorted       //
//   (2) size - Size of the input array                      //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Skewness of the array values. NaN will be returned for  //
//   NaN-only input arrays.                                  //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the skewness of the input data array values. //
//   Skewness is defined as                                  //
//                                                           //
//     mom3 / sqrt(mom2^3)                                   //
//                                                           //
//   where mom2 and mom3 are the second and third moment of  //
//   the data as returned by the moment() and moments()      //
//   functions. The skewness of normally distributed data    //
//   with a mean of zero should be 0.                        //
// --------------------------------------------------------- //

double skewness_flt(const float *data, const size_t size)
{
	double m2, m3, m4;
	moments_flt(data, size, mean_flt(data, size), &m2, &m3, &m4);
	return m3 / sqrt(m2 * m2 * m2);
}



// --------------------------------------------------------- //
// Kurtosis                                                  //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the data array to be sorted       //
//   (2) size - Size of the input array                      //
//                                                           //
// Returns:                                                  //
//                                                           //
//   Kurtosis of the array values. NaN will be returned for  //
//   NaN-only input arrays.                                  //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the kurtosis of the input data array values. //
//   Kurtosis is defined as                                  //
//                                                           //
//     mom4 / mom2^2                                         //
//                                                           //
//   where mom2 and mom4 are the second and fourth moment of //
//   the data as returned by the moment() and moments()      //
//   functions. The kurtosis of normally distributed data    //
//   with a mean of zero should be 3.                        //
// --------------------------------------------------------- //

double kurtosis_flt(const float *data, const size_t size)
{
	double m2, m3, m4;
	moments_flt(data, size, mean_flt(data, size), &m2, &m3, &m4);
	return m4 / (m2 * m2);
}



// --------------------------------------------------------- //
// Skewness and kurtosis                                     //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) data - Pointer to the data array to be sorted       //
//   (2) size - Size of the input array                      //
//   (3) skew - Pointer to variable for returning skewness   //
//   (4) kurt - Pointer to variable for returning kurtosis   //
//                                                           //
// Returns:                                                  //
//                                                           //
//   No return value.                                        //
//                                                           //
// Description:                                              //
//                                                           //
//   Calculates the skewness and kurtosis of the input data  //
//   array simultaneously. This is faster than calculating   //
//   them individually in situations where both are needed.  //
//   See the skewness() and kurtosis() functions for more    //
//   details on how these parameters are defined. Both skew  //
//   and kurt will be set to NaN if the array only contains  //
//   NaN values.                                             //
// --------------------------------------------------------- //

void skew_kurt_flt(const float *data, const size_t size, double *skew, double *kurt)
{
	double m2, m3, m4;
	moments_flt(data, size, mean_flt(data, size), &m2, &m3, &m4);
	*skew = m3 / sqrt(m2 * m2 * m2);
	*kurt = m4 / (m2 * m2);
	return;
}



// --------------------------------------------------------- //
// 1D boxcar filter                                          //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1)          data - Pointer to data array to be         //
//                       filtered.                           //
//   (2)     data_copy - Pointer to data array to be used    //
//                       for storing a copy of the data dur- //
//                       ring filtering. Its size must be    //
//                       equal to size + 2 * filter_radius.  //
//   (3)          size - Size of input array.                //
//   (4) filter_radius - Radius of boxcar filter.            //
//                                                           //
// Returns:                                                  //
//                                                           //
//   No return value.                                        //
//                                                           //
// Description:                                              //
//                                                           //
//   Applies a boxcar filter to the data array. NOTE that    //
//   this will modify the original data array. NaN values    //
//   will be set to 0 prior to filtering only if replace_nan //
//   is set to true (this is to avoid a check within this    //
//   function for reasons of speed). Values outside of the   //
//   boundaries of the array are assumed to be 0.            //
// --------------------------------------------------------- //

void filter_boxcar_1d_flt(float *data, float *data_copy, const size_t size, const size_t filter_radius)
{
	// Define filter size
	const size_t filter_size = 2 * filter_radius + 1;
	const float inv_filter_size = 1.0 / filter_size;
	size_t i;
	
	// Make copy of data, taking care of NaN
	for(i = size; i--;) data_copy[filter_radius + i] = FILTER_NAN(data[i]);

	// Fill overlap regions with 0
	for(i = filter_radius; i--;) data_copy[i] = data_copy[size + filter_radius + i] = 0.0;
	
	// Apply boxcar filter to last data point
	data[size - 1] = 0.0;
	for(i = filter_size; i--;) data[size - 1] += data_copy[size + i - 1];
	data[size - 1] *= inv_filter_size;
	
	// Recursively apply boxcar filter to  all previous data points
	for(i = size - 1; i--;) data[i] = data[i + 1] + (data_copy[i] - data_copy[filter_size + i]) * inv_filter_size;

	//int k=0;
	//data[size-1]=0.0;
	return;
}



// --------------------------------------------------------- //
// 2D Gaussian filter                                        //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1)          data - Pointer to data array to be         //
//                       filtered.                           //
//   (2)     data_copy - Pointer to data array to be used    //
//                       for storing a single column of the  //
//                       data array. Its size must be equal  //
//                       to size_y.                          //
//   (3)      data_row - Pointer to data array to be used by //
//                       the boxcar filter to be employed.   //
//                       Its size must be equal to size_x +  //
//                       2 * filter_radius.                  //
//   (4)      data_col - Pointer to data array to be used by //
//                       the boxcar filter to be employed.   //
//                       Its size must be equal to size_y +  //
//                       2 * filter_radius.                  //
//   (5)        size_x - Size of the first dimension of the  //
//                       input data array.                   //
//   (6)        size_y - Size of the second dimension of the //
//                       input data array.                   //
//   (7)        n_iter - Number of iterations of boxcar fil- //
//                       tering to be carried out.           //
//   (8) filter_radius - Radius of the boxcar filter to be   //
//                       applied. The filter width will be   //
//                       defined as 2 * filter_radius + 1.   //
//                                                           //
// Returns:                                                  //
//                                                           //
//   No return value.                                        //
//                                                           //
// Description:                                              //
//                                                           //
//   Applies a pseudo-Gaussian filter to the two-dimensional //
//   data array by running a series of n_iter boxcar filters //
//   of radius filter_radius across the data array in both   //
//   dimensions. The function optimal_filter_size(...) can   //
//   be used to automatically determine the required number  //
//   of iterations and boxcar filter radius for a given      //
//   standard deviation of the Gaussian.                     //
//                                                           //
//   For reasons of speed several data arrays must have been //
//   pre-allocated and passed on to this function:           //
//                                                           //
//   - data_copy: Used to store a single column of the input //
//                data. Must be of size size_y.              //
//   - data_row:  Used to store a copy of the data passed on //
//                to the boxcar filter. Must be of size      //
//                size_x + 2 * filter_radius.                //
//   - data_col:  Used to store a copy of the data passed on //
//                to the boxcar filter. Must be of size      //
//                size_y + 2 * filter_radius.                //
//                                                           //
//   The sole purpose of having these array created extern-  //
//   ally and then passed on to the function is to improve   //
//   the speed of the algorithm in cases where it needs to   //
//   be invoked a large number of times.                     //
//                                                           //
//   Note that the function will only be able to approximate //
//   a Gaussian, with the approximation becoming better for  //
//   a larger number of iterations of the boxcar filter. 3-4 //
//   iterations should already provide a reasonable approxi- //
//   mation for the purpose of image filtering. Also note    //
//   that the value of the standard deviation of the Gauss-  //
//   ian will only be approximated by the optimal number of  //
//   iterations and filter radius, usually to within several //
//   percent.                                                //
// --------------------------------------------------------- //

void filter_gauss_2d_flt(float *data, float *data_copy, float *data_row, float *data_col, const size_t size_x, const size_t size_y, const size_t n_iter, const size_t filter_radius)
{
	// Set up a few variables
	const size_t size_xy = size_x * size_y;
	float *ptr = data + size_xy;
	float *ptr2;
	float *ptr_cpy=ptr;


	// Run row filter (along x-axis)
	// This is straightforward, as the data are contiguous in x.
	while(ptr > data)
	{
		ptr -= size_x;
		for(size_t i = n_iter; i--;){
			 filter_boxcar_1d_flt(ptr, data_row, size_x, filter_radius);
		}
	}


	// Run column filter (along y-axis)
	// This is more complicated, as the data are non-contiguous in y.
	for(size_t x = size_x; x--;)
	{
		// Copy data into column array
		ptr = data + size_xy - size_x + x;
		ptr2 = data_copy + size_y;
		while(ptr2 --> data_copy)
		{
			*ptr2 = *ptr;
			ptr -= size_x;
		}
		
		// Apply filter
		for(size_t i = n_iter; i--;) filter_boxcar_1d_flt(data_copy, data_col, size_y, filter_radius);
		
		// Copy column array back into data array
		ptr = data + size_xy - size_x + x;
		ptr2 = data_copy + size_y;
		while(ptr2 --> data_copy)
		{
			*ptr = *ptr2;
			ptr -= size_x;
		}
	}
	return;
}



// --------------------------------------------------------- //
// Determine optimal boxcar filter size for Gaussian filter  //
// --------------------------------------------------------- //
//                                                           //
// Arguments:                                                //
//                                                           //
//   (1) sigma         - Standard deviation of the Gaussian. //
//   (2) filter_radius - Radius of the boxcar filter that is //
//                       to be used in approximation.        //
//   (3) n_iter        - Number of iteration of boxcar fil-  //
//                       tering to be used in approximation. //
//                                                           //
// Returns:                                                  //
//                                                           //
//   No return value.                                        //
//                                                           //
// Description:                                              //
//                                                           //
//   Based on the specified standard deviation of the Gauss- //
//   ian filter, this function will determine the optimal    //
//   radius and number of iterations for the boxcar filter   //
//   that is to be used to approximate the Gaussian kernel.  //
//   Both 'filter_radius' and 'n_iter' are pointers and will //
//   be set to the correct values by this function. They can //
//   then be fed into the Gaussian filtering function.       //
// --------------------------------------------------------- //

void optimal_filter_size_flt(const double sigma, size_t *filter_radius, size_t *n_iter)
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
