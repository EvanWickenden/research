
#include <random>
#include <iostream>
#include <cstdio>
#include <stdexcept>

#include "su2.h"
#include "lattice.h"




void Lattice::store(FILE *fp)
{
    if (fwrite(dimensions, sizeof (int), 4, fp) < 4)
    {
        std::cerr << "fwrite() failed" << std::endl;
        exit(1);
    }
    int size = t_increment * dimensions[0];
    if (fwrite(sites, sizeof (LatticeSite), size, fp) < size)
    {
        std::cerr << "fwrite() failed" << std::endl;
        exit(1);
    }
}

void Lattice::load(FILE *fp)
{
    if (fread(dimensions, sizeof (int), 4, fp))
    {
        std::cerr << "fread() failed" << std::endl;
        exit(1);
    }
    t_increment = dimensions[1] * dimensions[2] * dimensions[3];
    x_increment = dimensions[2] * dimensions[3];
    y_increment = dimensions[3];
     
    acceptance_ratio[0] = 0;
    acceptance_ratio[1] = 0;

    int size = t_increment * dimensions[0]; 
    if (fread(sites, sizeof (LatticeSite), size, fp) < size)
    {
        std::cerr << "fread() failed" << std::endl;
        exit(1);
    }
}

std::ostream& operator<<(std::ostream& stream, Lattice::LatticeSite& site)
{
    stream << "lattice site: " << std::endl;
    stream << "forward_t:\n" << site.forward_t << std::endl;
    stream << "forward_x:\n" << site.forward_x << std::endl;
    stream << "forward_y:\n" << site.forward_y << std::endl;
    stream << "forward_z:\n" << site.forward_z << std::endl;
    return stream;
}

void Lattice::re_unitarize()
{
    int i;
    for (i = 0; i < t_increment * dimensions[0]; i++)
    {
        sites[i].forward_t.re_unitarize();
        sites[i].forward_x.re_unitarize();
        sites[i].forward_y.re_unitarize();
        sites[i].forward_z.re_unitarize();
    }
}


Lattice::Lattice(const int dimensions[4], LatticeSite *sites) :
    t_increment(dimensions[1] * dimensions[2] * dimensions[3]),
    x_increment(dimensions[2] * dimensions[3]),
    y_increment(dimensions[3]),
    sites(new LatticeSite[t_increment * dimensions[0]])
{
    this->dimensions[0] = dimensions[0];
    this->dimensions[1] = dimensions[1];
    this->dimensions[2] = dimensions[2];
    this->dimensions[3] = dimensions[3];

    acceptance_ratio[0] = 0;
    acceptance_ratio[1] = 0;

    std::memcpy(this->sites, sites, t_increment * dimensions[0] * sizeof (LatticeSite));
}

Lattice::Lattice(const int dimensions[4], std::mt19937& gen, double width) : 
    t_increment(dimensions[1] * dimensions[2] * dimensions[3]),
    x_increment(dimensions[2] * dimensions[3]),
    y_increment(dimensions[3]),
    sites(new LatticeSite[dimensions[0] * dimensions[1] * dimensions[2] * dimensions[3]])
{
    this->dimensions[0] = dimensions[0];
    this->dimensions[1] = dimensions[1];
    this->dimensions[2] = dimensions[2];
    this->dimensions[3] = dimensions[3];

    acceptance_ratio[0] = 0;
    acceptance_ratio[1] = 0;

    int i;
    for (i = 0; i < t_increment * dimensions[0]; i++)
    {
        sites[i].forward_t = su2::generate2(gen, width); 
        sites[i].forward_x = su2::generate2(gen, width); 
        sites[i].forward_y = su2::generate2(gen, width); 
        sites[i].forward_z = su2::generate2(gen, width); 
    }
}

Lattice::~Lattice()
{
    delete[] sites;
}


Lattice::LatticeSite* Lattice::operator()(int t, int x, int y, int z)
{
    t = (t + dimensions[0]) % dimensions[0];
    x = (x + dimensions[1]) % dimensions[1];
    y = (y + dimensions[2]) % dimensions[2];
    z = (z + dimensions[3]) % dimensions[3];

    return &sites[t*t_increment + x*x_increment + y*y_increment + z];
}



double Lattice::total_action()
{
    int i,j,k,l;

    double sum = 0;

    for (i = 0; i < dimensions[0]; i++)
    {
     for (j = 0; j < dimensions[0]; j++)
     {
      for (k = 0; k < dimensions[0]; k++)
      {
       for (l = 0; l < dimensions[0]; l++)
       {
           sum += tx_plaquette(i,j,k,l).trace().real()
               + ty_plaquette(i,j,k,l).trace().real()
               + tz_plaquette(i,j,k,l).trace().real()
               + xy_plaquette(i,j,k,l).trace().real()
               + xz_plaquette(i,j,k,l).trace().real()
               + yz_plaquette(i,j,k,l).trace().real();
       }
      }
     }
    }
    return sum;
}


su2 Lattice::tx_plaquette(int t, int x, int y, int z)
{
    LatticeSite *bottom_left((*this)(t, x, y, z));
    LatticeSite *top_left((*this)(t + 1, x, y, z));
    LatticeSite *bottom_right((*this)(t, x + 1, y, z));

    return bottom_left->forward_t.inverse()
        * top_left->forward_x.inverse()
        * bottom_right->forward_t
        * bottom_left->forward_x;
}

su2 Lattice::ty_plaquette(int t, int x, int y, int z)
{
    LatticeSite *bottom_left((*this)(t, x, y, z));
    LatticeSite *top_left((*this)(t + 1, x, y, z));
    LatticeSite *bottom_right((*this)(t, x, y + 1, z));

    return bottom_left->forward_t.inverse()
        * top_left->forward_y.inverse()
        * bottom_right->forward_t
        * bottom_left->forward_y;
}

su2 Lattice::tz_plaquette(int t, int x, int y, int z)
{
    LatticeSite *bottom_left((*this)(t, x, y, z));
    LatticeSite *top_left((*this)(t + 1, x, y, z));
    LatticeSite *bottom_right((*this)(t, x, y, z + 1));

    return bottom_left->forward_t.inverse()
        * top_left->forward_z.inverse()
        * bottom_right->forward_t
        * bottom_left->forward_z;
}

su2 Lattice::xy_plaquette(int t, int x, int y, int z)
{
    LatticeSite *bottom_left((*this)(t, x, y, z));
    LatticeSite *top_left((*this)(t, x + 1, y, z));
    LatticeSite *bottom_right((*this)(t, x, y + 1, z));

    return bottom_left->forward_x.inverse()
        * top_left->forward_y.inverse()
        * bottom_right->forward_x
        * bottom_left->forward_y;
}

su2 Lattice::xz_plaquette(int t, int x, int y, int z)
{
    LatticeSite *bottom_left((*this)(t, x, y, z));
    LatticeSite *top_left((*this)(t, x + 1, y, z));
    LatticeSite *bottom_right((*this)(t, x, y, z + 1));

    return bottom_left->forward_x.inverse()
        * top_left->forward_z.inverse()
        * bottom_right->forward_x
        * bottom_left->forward_z;
}

su2 Lattice::yz_plaquette(int t, int x, int y, int z)
{
    LatticeSite *bottom_left((*this)(t, x, y, z));
    LatticeSite *top_left((*this)(t, x, y + 1, z));
    LatticeSite *bottom_right((*this)(t, x, y, z + 1));

    return bottom_left->forward_y.inverse()
        * top_left->forward_z.inverse()
        * bottom_right->forward_y
        * bottom_left->forward_z;
}

double Lattice::t_link_action(int t, int x, int y, int z)
{
    double sum = 0;
    sum += tx_plaquette(t, x, y, z).trace().real()
        + tx_plaquette(t, x - 1, y, z).trace().real()
        + ty_plaquette(t, x, y, z).trace().real()
        + ty_plaquette(t, x, y - 1, z).trace().real()
        + tz_plaquette(t, x, y, z).trace().real()
        + tz_plaquette(t, x, y, z - 1).trace().real();
    return sum;
}

double Lattice::x_link_action(int t, int x, int y, int z)
{
    double sum = 0;
    sum += tx_plaquette(t, x, y, z).trace().real()
        + tx_plaquette(t - 1, x, y, z).trace().real()
        + xy_plaquette(t, x, y, z).trace().real()
        + xy_plaquette(t, x, y - 1, z).trace().real()
        + xz_plaquette(t, x, y, z).trace().real()
        + xz_plaquette(t, x, y, z - 1).trace().real();
    return sum;
}

double Lattice::y_link_action(int t, int x, int y, int z)
{
    double sum = 0;
    sum += ty_plaquette(t, x, y, z).trace().real()
        + ty_plaquette(t - 1, x, y, z).trace().real()
        + xy_plaquette(t, x, y, z).trace().real()
        + xy_plaquette(t, x - 1, y, z).trace().real()
        + yz_plaquette(t, x, y, z).trace().real()
        + yz_plaquette(t, x, y, z - 1).trace().real();
    return sum;
}

double Lattice::z_link_action(int t, int x, int y, int z)
{
    double sum = 0;
    sum += tz_plaquette(t, x, y, z).trace().real()
        + tz_plaquette(t - 1, x, y, z).trace().real()
        + xz_plaquette(t, x, y, z).trace().real()
        + xz_plaquette(t, x - 1, y, z).trace().real()
        + yz_plaquette(t, x, y, z).trace().real()
        + yz_plaquette(t, x, y - 1, z).trace().real();
    return sum;
}

