/*
 * Project: 3D Total Variation minimization
 * Author: Gokhan Gunay, ghngunay@gmail.com
 * Copyright: (C) 2018 by Gokhan Gunay
 * License: GNU GPL v3 (see License.txt)
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <new>
#include <math.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "ArgumentParser.hpp"
#include "tv_filter.h"

using namespace std;

int main(int argc, char * argv[])
{
	Cparser parser(argc, argv);
	parser.save_key("in_file", "-in");
	parser.save_key("out_file", "-out");
	parser.save_key("lambda", "-l");
	parser.save_key("iter", "-it");
	parser.save_key("verbose", "-v");
	parser.save_key("IsIsotropic", "-iso");
	parser.save_key("SliceBySlice", "-slc");

	typedef itk::Image<float, 3> InputImageType;
	typedef itk::ImageFileReader<InputImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	typedef itk::TotalVariationMinimization<InputImageType, InputImageType> TV;
	TV::Pointer Tv = TV::New();

	if (parser["in_file"].is_called())
	{
		reader->SetFileName(parser["in_file"].get_as_string()[0]);		
		cout << "In File:" << parser["in_file"].get_as_string()[0] << "\n";
		Tv->SetInput(reader->GetOutput());
		Tv->SetIsotropic(parser["IsIsotropic"].is_called());
		Tv->SetSliceBySlice(parser["SliceBySlice"].is_called());

		auto lmbd = parser["lambda"].get_as_float();
		if (!std::empty(lmbd))
		{
			Tv->SetLambda(lmbd[0]);
		}
		auto it = parser["iter"].get_as_integer();
		if (!std::empty(it))
		{
			Tv->SetIt(it[0]);
		}
	}
	else
	{
		std::cout << "No input file set!!" << std::endl;
	}
	if (parser["out_file"].is_called())
	{
		typedef itk::ImageFileWriter<InputImageType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(Tv->GetOutput());
		writer->SetFileName(parser["out_file"].get_as_string()[0]);
		cout << "Out File:" << parser["out_file"].get_as_string()[0] << "\n";
		try
		{
			writer->Update();
		}
		catch (...)
		{
			std::cerr << "Invalid output image" << std::endl;
			return EXIT_FAILURE;
		}
	}
	else
	{
		std::cout << "No output file set!!" << std::endl;
	}

	return 0;
}
