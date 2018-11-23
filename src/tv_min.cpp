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


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"



#include "./Parser/ArgumentParser.hpp"
#include <math.h>

#include "tv.h"

using namespace std;


int main(int argc,char * argv[])
{
	float lambda=0;
	unsigned int It=20;
	bool verbose=false;
	cparser parser(argc,argv);
	parser.save_keys("in_file", "-in");
	parser.save_keys("out_file", "-out");
	parser.save_keys("lambda", "-l");
	parser.save_keys("iter", "-it");
	parser.save_keys("verbose", "-v");
	parser.save_keys("slicebyslice", "-slc");
	cout<<"In File:"<<parser["in_file"]<<"\n";
	cout<<"In File:"<<parser["out_file"]<<"\n";
	try
	{
		lambda=stof(parser["lambda"]);
	}
	catch(const std::invalid_argument& ia)
    {
    	std::cout << "Invalid lambda value 0 is assigned instead\n";
    }
    try
	{
		It=stof(parser["iter"]);
	}
	catch(const std::invalid_argument& ia)
    {
    	std::cout << "Invalid iteration value 10 is assigned instead\n";
    }
	

	typedef itk::Image< float, 3> _InputImageType;
	typedef itk::ImageFileReader< _InputImageType > _ReaderType;
	_ReaderType::Pointer reader = _ReaderType::New();

	reader->SetFileName( parser["in_file"] );
    
    typedef itk::TotalVariationMinimization<_InputImageType, _InputImageType> TV_;

	TV_::Pointer TV=TV_::New();
	TV->SetInput(reader->GetOutput());
	if (parser["slicebyslice"] == "false")
		TV->SlcBySlc(false);
	TV->SetIt(It);
	TV->SetLambda(lambda);

	typedef itk::ImageFileWriter< _InputImageType > _WriterType;
	_WriterType::Pointer writer = _WriterType::New();
	writer->SetInput(TV->GetOutput());
	writer->SetFileName( parser["out_file"] );
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & err )
    {
    	std::cerr << "Invalid output image" << std::endl;
    	return EXIT_FAILURE;
    }
}
