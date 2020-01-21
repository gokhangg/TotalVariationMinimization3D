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

#include "../CommandlineParser/ArgumentParser.hpp"
#include "tv_filter.h"

using namespace std;

int main(int argc, char * argv[])
{
    float lambda = 0;
    unsigned int it = 20;
    bool verbose = false;
    cparser parser(argc, argv);
    parser.save_keys("in_file", "-in");
    parser.save_keys("out_file", "-out");
    parser.save_keys("lambda", "-l");
    parser.save_keys("iter", "-it");
    parser.save_keys("verbose", "-v");
    parser.save_keys("IsIsotropic", "-iso");
	parser.save_keys("SliceBySLice", "-slc");
    cout << "In File:" << parser["in_file"] << "\n";
    cout << "In File:" << parser["out_file"] << "\n";
    try
    {
        lambda = stof(parser["lambda"]);
    } 
    catch (const std::invalid_argument& ia)
    {
        std::cout << "Invalid lambda value 0 is assigned instead\n";
    }
    try
    {
        it = stof(parser["iter"]);
    } 
    catch (const std::invalid_argument& ia)
    {
        std::cout << "Invalid iteration value 20 is assigned instead\n";
    }

    typedef itk::Image<float, 3> InputImageType;
    typedef itk::ImageFileReader<InputImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(parser["in_file"]);

    typedef itk::TotalVariationMinimization<InputImageType, InputImageType> TV;

    TV::Pointer Tv = TV::New();
    Tv->SetInput(reader->GetOutput());
    if (parser["IsIsotropic"] == "true")
    {
        Tv->SetIsotropic(true);
    }
	if (parser["SliceBySlice"] == "true")
	{
		Tv->SetSliceBySlice(true);
	}
    Tv->SetIt(it);
    Tv->SetLambda(lambda);

    typedef itk::ImageFileWriter<InputImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(Tv->GetOutput());
    writer->SetFileName(parser["out_file"]);
    try
    {
        writer->Update();
    } catch (itk::ExceptionObject & err)
    {
        std::cerr << "Invalid output image" << std::endl;
        return EXIT_FAILURE;
    }

	return 0;
}
