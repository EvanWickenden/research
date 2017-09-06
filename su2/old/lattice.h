#ifndef __LATTICE_H__
#define __LATTICE_H__

#include <iostream>

#include "su2.h"


struct Lattice
{
    struct LatticeSite
    {
        su2 forward_t;
        su2 forward_x;
        su2 forward_y;
        su2 forward_z;

        friend std::ostream& operator<<(std::ostream& stream, LatticeSite& site);
    };

    int dimensions[4]; /* t, x, y, z */

    int t_increment;
    int x_increment;
    int y_increment;

    int acceptance_ratio[2];

    LatticeSite *sites;


    Lattice(const int dimensions[4], LatticeSite *sites);
    Lattice(const int dimensions[4], std::mt19937& gen, double width = 0.5); 

    ~Lattice();

    LatticeSite *operator()(int t, int x, int y, int z);

    void re_unitarize();

    void store(FILE *fp);
    void load(FILE *fp);

    double total_action();
    
    su2 tx_plaquette(int t, int x, int y, int z);
    su2 ty_plaquette(int t, int x, int y, int z);
    su2 tz_plaquette(int t, int x, int y, int z);
    su2 xy_plaquette(int t, int x, int y, int z);
    su2 xz_plaquette(int t, int x, int y, int z);
    su2 yz_plaquette(int t, int x, int y, int z);

    double t_link_action(int t, int x, int y, int z);
    double x_link_action(int t, int x, int y, int z);
    double y_link_action(int t, int x, int y, int z);
    double z_link_action(int t, int x, int y, int z);

};



#endif 
