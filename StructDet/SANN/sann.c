/*
* @file   sann.c
* @author van Meel, Filion, Valeriani, Frenkel
* @date   November 2011
* @brief  Sample implementation of the SANN algorithm in C
*/


/*
* @struct NbData
* @brief Defines an {index,distance} pair for use as neighbour data
*/
struct NbData
{
	int   id; // Id of neighbour particle
	double distance; // Distance to neighbour particle
};

/*
* @brief Compares the distance of two neighbours for use in 'qsort'
* @param nb1 Pointer to first neighbour's NbData
* @param nb2 Pointer to second neighbour's NbData
* @retval Returns negative if less, positive otherwise
*/
int nbLess( const void *nb1, const void *nb2 )
{
	// cast 'const void' pointer to 'const NbData' pointer for both neighbours
	NbData *pnb1 = (const NbData*) nb1;
	NbData *pnb2 = (const NbData*) nb2;

	// If first distance is smaller than second distance return negative value
	if ( pnb1->distance < pnb2->distance ) return -1;

	// Forget about 'equal' when using floating point numbers...
	// Return positive value for 'greater than'
	return 1;
}

/*
* @brief Computes the SANN set of nearest neighbours for a given particle
* @param id Id of particle who's neighbours are to be computed
* @param neighbours NbData array to receive neighbour {id,distance} pairs
* @param radius Pointer to double to receive SANN radius
* @retval Number of neighbours computed by SANN, (-1) on error
*/
int computeSANN ( int id, NbData *neighbours, double *radius )
{
	double distanceSum; // sum of neighbour distances
	int count; // number of all potential neighbours available
	int i; // a loop variable

	//Step 1:
	//   Get number and {id,distance} pairs of all potential neighbours.
	//   In this example we use a Verlet neighbour list with a
	//   large-enough cutoff distance for this task. SANN then chooses
	//   its neighbours from this set.
	count = computeVerletNeighbors( id, neighbours );

	// If there are not enough neighbours available, report an error
	if (count < 3) return -1;

	// Step 2:
	//   Sort neighbours according to their distance in increasing order.
	//   In this example we use a 'quicksort' algorithm for this task,
	//   which exists in the standard C library stdlib.h
	qsort( neighbours, count, sizeof( NbData ), nbLess );

	// Step 3 / 4:
	//   Start with 3 neighbours (it's the minimum number of
	//   neighbours possible)
	distanceSum = 0;
	for (i=0; i<3; ++i)
	{
		// Add neighbour distance to sum
		distanceSum += nbData[i].distance;
	}
	// Set SANN radius to distanceSum / (i - 2)
	*radius = distanceSum;

	// Step 4 / 5:
	//   Iteratively include further neighbours until finished, which is if
	//   the SANN radius is smaller than the distance to the next neighbour
	while ((i < count) && (radius > neighbours[i].distance))
	{
		// Add neighbour distance to sum
		distanceSum += neighbours[i].distance;
		// Compute new SANN radius
		*radius = distanceSum / (i - 2.0);
		// increase the SANN number of neighbours
		++i;
	}

	// If there were not enough neighbours for the algorithm to converge,
	// report an error
	if (i == count) return -1;

	// Step 6:
	//   Return the number of SANN neighbours.
	//   Note: the SANN radius has already been stored in the pointer 'radius',
	//   which was provided as parameter to the function
	return i;
}

// end-of-file
