#ifndef TESTS_H
#define TESTS_H

/* A class for testing Brute Force optimization method */

namespace BioLCCC {

class BruteForceTester {
	public:
                BruteForceTester ();
		double calculate ();

		void set_x (double new_x);
		void set_y (double new_y);
                void set_z (double new_z);

	private:
                double x, y, z;
};

}

#endif
