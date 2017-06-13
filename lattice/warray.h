#ifndef __WARRAY_H__
#define __WARRAY_H__

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "log.h"

template <typename T>
struct WArray
{
	T* data;
	int N;

	WArray(int N)
	{
		data = (T *) malloc(N * sizeof (T));
		if (data == NULL)
			DIE("malloc() returned NULL");
	}

	~WArray()
	{
		free(data);
	}

	T& operator[](int i) { return data[i]; }
	const T& operator[](int i) const { return data[i]; }

	WArray(const char *load_filename)
	{
		int objects_read;
		FILE *fp = fopen(load_filename, "r");
		if (fp == NULL)
			DIE("fopen() failed");

		objects_read = fread(&N, sizeof (int), 1, fp);
		if (objects_read != 1)
			DIE("fread() failed");

		data = (T *) malloc(N * sizeof (T));
		if (data == NULL)
			DIE("malloc() returned NULL");

		objects_read = fread(data, sizeof (T), N, fp);
		if (objects_read != N)
			DIE("fread() failed");
	}

	void record(const char *save_filename)
	{
		int objects_written;
		FILE *fp = fopen(save_filename, "w");
		if (fp == NULL)
			DIE("fopen() failed");

		objects_written = fwrite(&N, sizeof (int), 1, fp);
		if (objects_written != 1)
			DIE("fwrite() failed");

		objects_written = fwrite(data, sizeof (T), N, fp);
		if (objects_written != N)
			DIE("fwrite() falide");
	}


	friend std::ostream& operator<<(std::ostream& s, const WArray<T>& arr)
	{
		int i;
		s << "[";
		for (i = 0; i < arr.N - 1; i++)
		{
			s << arr[i] << ", ";
		}
		s << arr[i] << "]" << std::endl;
		return s;
	}
};


#endif
