/*
 * Project: 3D Total Variation minimization
 * Author: Gokhan Gunay, ghngunay@gmail.com
 * Copyright: (C) 2018 by Gokhan Gunay
 * License: GNU GPL v3 (see License.txt)
 */

#include "tv.h"

#ifndef tv_hxx
#define tv_hxx

#define epsilon 0.0000001

namespace itk
{

	template<typename TInputImage, typename TOutputImage>
	void
	TotalVariationMinimization<TInputImage, TOutputImage>
	::GenerateData()
	{
		//Allocate
		uint32 Dim=TInputImage::ImageDimension;
		float *pP[3];
		float *pIm;
		float r = 0.0;
		typename TInputImage::Pointer InIm = const_cast< TInputImage * >(this->GetInput());

		/*If the image is 3D and desired to be processed slice by slice in X-Y plane*/
		if (SliceBySlice)
		Dm = 2;
		/*If the image is 3D and desired to be processed as 3D*/
		else
		Dm = Dim;

		for (int i = 0; i<Dim; i++)
		{
			Sz[i] = InIm->GetBufferedRegion().GetSize()[i];
		}

		for (int i = 0; i<Dm; i++)
		{
			TotalVox *= Sz[i] + 2;
			Sp[i] = InIm->GetSpacing()[i];
		}

		if (!IsIsotropic)
		{
			for (int i = 0; i < Dm; i++)
			r += 1 / ((Sp[i] + epsilon)*(Sp[i] + epsilon));
			r = sqrt(Dm / r);
			for (int i = 0; i < Dm; i++)
			Sc[i] = r / (Sp[i] + epsilon);
		}
		else
		{
			for (int i = 0; i < Dim; i++)
			Sc[i] = 1;
		}

		for (int i = 0; i<Dm; i++)
		{
			P[i] = new InPxType[TotalVox];
			pP[i] = P[i];
		}

		Im = new InPxType[TotalVox];
		Allocated = true;

		pIm = Im;
		for (int i = (TotalVox); i >= 0; i--)
		{
			for (int j = 0; j<Dm; j++)
			*(pP[j]++) = 0;
			*(pIm++) = 0;
		}
		ImLoaded = true;

		this->GetOutput()->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
		this->GetOutput()->Allocate();

		{ //Update
			if (Dim == 2)
			{
				this->EngineSliceBySlice();
			}
			else if (Dim == 3)
			{
				if (SliceBySlice)
				{
					this->EngineSliceBySlice();
				}
				else
				{
					this->Load3D();
					this->Engine3D();
				}
			}
		}
	}

	template<typename TInputImage, typename TOutputImage>
	void
	TotalVariationMinimization<TInputImage, TOutputImage>
	::Load3D(void)
	{
		uint32 Dim =TInputImage::ImageDimension;
		InPxType *pIm;
		uint32 Sz2[3];
		uint32 SzPl;
		ImageRegionConstIterator<TInputImage> InIt(const_cast< TInputImage *>(this->GetInput()),
				const_cast<TInputImage*>(this->GetInput())->GetRequestedRegion());

		for(int i=0;i<Dm;i++)
		Sz2[i]=Sz[i]+2;
		SzPl = Sz2[0] * Sz2[1];

		InIt.GoToBegin();
		if (Dm == 3)
		{
			pIm = Im + SzPl + Sz2[0] + 1;
			for (int k = 0; k < Sz[2]; k++)
			{
				for (int j = 0; j < Sz[1]; j++)
				{
					for (int i = 0; i < Sz[0]; i++)
					{
						*pIm++ = static_cast<InPxType>(InIt.Get());
						++InIt;
					}
					pIm += 2;
				}
				pIm += 2 * Sz2[0];
			}
		}
	}

	template<typename TInputImage, typename TOutputImage>
	void
	TotalVariationMinimization<TInputImage, TOutputImage>
	::Load2D(void)
	{
		uint32 Dim =TInputImage::ImageDimension;
		InPxType *pIm;
		uint32 Sz2[3];
		typename TInputImage::Pointer InIm = const_cast< TInputImage * >(this->GetInput());
		ImageRegionConstIterator<TInputImage> InIt(InIm, InIm->GetRequestedRegion());
		for(int i=0;i<Dm;i++)
		Sz2[i]=Sz[i]+2;
		InIt.GoToBegin();
		if (Dm == 2)
		{
			pIm = Im + Sz2[0] + 1;
			for (int j = 0; j < Sz[1]; j++)
			{
				for (int i = 0; i < Sz[0]; i++)
				{
					*pIm++ = static_cast<InPxType>(InIt.Get());
					++InIt;
				}
				pIm += 2;
			}
		}
	}

//Reads from In image, applies 2D TV minimization and saves the result to Out back.
//In order to reduce memory use, all gradient and divergence operators applied together in the for loops.
	template<typename TInputImage,typename TOutputImage>
	void
	TotalVariationMinimization<TInputImage, TOutputImage>
	::TVmin2D()
	{
		uint32 Dim=TInputImage::ImageDimension;
		uint32 Sz2[3];
		float Lam = Lm + epsilon;
		float rat;
		InPxType *pP[3];
		InPxType *pIm;
		InPxType diff[3];
		typename TOutputImage::Pointer OutIm = this->GetOutput();

		if(ImLoaded)
		{
			for(int i=0;i<Dm;i++)
			Sz2[i]=Sz[i]+2;
			for (int it=0;it<It;it++)
			{
				for(int i=0;i<Dm;i++)
				{
					pP[i]=P[i]+Sz2[0]+1;
				}
				pIm=Im+Sz2[0]+1;
				for(int j=0;j<Sz[1];j++)
				{
					for(int i=0;i<Sz[0];i++)
					{
						diff[0] = Sc[0] * ((*pP[0]) - *(pP[0] - 1)) + Sc[1] * (*(pP[1]) - *(pP[1] - Sz2[0])) - *(pIm) / Lam;
						diff[1]=diff[0];
						diff[0] = Sc[0] * (Sc[0] * (*(pP[0] + 1) - *(pP[0])) + Sc[1] * (*(pP[1] + 1) - *(pP[1] - Sz2[0] + 1)) - *(pIm + 1) / Lam - diff[0]);
						diff[1] = Sc[1] * (Sc[0] * (*(pP[0] + Sz2[0]) - *(pP[0] - 1 + Sz2[0])) + Sc[1] * (*(pP[1] + Sz2[0]) - *(pP[1])) - *(pIm + Sz2[0]) / Lam - diff[1]);
						rat=(1+To*sqrt(diff[0]*diff[0]+diff[1]*diff[1]));
						*pP[0]=(*pP[0]+To*diff[0])/rat;
						*pP[1]=(*pP[1]+To*diff[1])/rat;
						for(int i=0;i<Dm;i++)
						pP[i]++;
						pIm++;
					}
					pIm+=2;
					for(int i=0;i<Dm;i++)
					pP[i]+=2;
				}
			}
			for(int i=0;i<Dm;i++)
			pP[i]=P[i]+1+(Sz2[0]);
			pIm=Im+1+(Sz2[0]);
			for(int j=0;j<Sz[1];j++)
			{
				for(int i=0;i<Sz[0];i++)
				{
					*pIm++ = (*pIm - Lam*(Sc[0] * (*pP[0] - *(pP[0] - 1)) + Sc[1] * (*pP[1] - *(pP[1] - Sz2[0]))));
					for(int k=0;k<Dm;k++)
					pP[k]++;
				}
				pIm+=2;
				for(int k=0;k<Dm;k++)
				pP[k]+=2;
			}
		}
		else
		{
		}
	}
//Reads from In image, applies 3D TV minimization and saves the result to Out back.
//In order to reduce memory use, all gradient and divergence operators applied together in the for loops.
	template<typename TInputImage, typename TOutputImage>
	void
	TotalVariationMinimization<TInputImage, TOutputImage>
	::Engine3D()
	{
		uint32 Dim=TInputImage::ImageDimension;
		uint32 Sz2[3];
		uint32 SzPl;
		float rat,Lam = Lm + epsilon;
		typename TOutputImage::Pointer OutIm = this->GetOutput();
		InPxType *pP[3];
		InPxType *pIm;
		InPxType diff[3];
		ImageRegionIterator<TOutputImage> OutIt(OutIm, OutIm->GetRequestedRegion());

		if (ImLoaded)
		{
			for (int i = 0; i<Dim; i++)
			Sz2[i] = Sz[i] + 2;
			SzPl = Sz2[0] * Sz2[1];

			for (int it = 0; it<It; it++)
			{
				for (int i = 0; i<Dim; i++)
				{
					pP[i] = P[i] + SzPl + Sz2[0] + 1;
				}
				pIm = Im + SzPl + Sz2[0] + 1;
				for (int k = 0; k < Sz[2]; k++)
				{
					for (int j = 0; j < Sz[1]; j++)
					{
						for (int i = 0; i < Sz[0]; i++)
						{
							diff[0] = Sc[0]*((*pP[0]) - *(pP[0] - 1)) + Sc[1]*(*(pP[1]) - *(pP[1] - Sz2[0])) + Sc[2]*(*(pP[2]) - *(pP[2] - SzPl)) - *(pIm) / Lam;
							diff[1] = diff[0];
							diff[2] = diff[0];
							diff[0] = Sc[0]*(Sc[0]*(*(pP[0] + 1) - *(pP[0])) + Sc[1]*(*(pP[1] + 1) - *(pP[1] - Sz2[0] + 1)) + Sc[2]*(*(pP[2] + 1) - *(pP[2] - SzPl + 1)) - *(pIm + 1) / Lam - diff[0]);
							diff[1] = Sc[1]*(Sc[0]*(*(pP[0] + Sz2[0]) - *(pP[0] - 1 + Sz2[0])) + Sc[1]*(*(pP[1] + Sz2[0]) - *(pP[1])) + Sc[2]*(*(pP[2] + Sz2[0]) - *(pP[2] - SzPl + Sz2[0])) - *(pIm + Sz2[0]) / Lam - diff[1]);
							diff[2] = Sc[2]*(Sc[0]*(*(pP[0] + SzPl) - *(pP[0] - 1 + SzPl)) + Sc[1]*(*(pP[1] + SzPl) - *(pP[1] + SzPl - Sz2[0])) + Sc[2]*(*(pP[2] + SzPl) - *(pP[2])) - *(pIm + SzPl) / Lam - diff[2]);
							rat = (1 + To*sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]));
							*pP[0] = (*pP[0] + To*diff[0]) / rat;
							*pP[1] = (*pP[1] + To*diff[1]) / rat;
							*pP[2] = (*pP[2] + To*diff[2]) / rat;
							for (int i = 0; i < Dim; i++)
							pP[i]++;
							pIm++;
						}
						pIm += 2;
						for (int i = 0; i < Dim; i++)
						pP[i] += 2;
					}
					pIm += 2*Sz2[0];
					for (int i = 0; i < Dim; i++)
					pP[i] += 2*Sz2[0];
				}
			}

			for (int l = 0; l<Dim; l++)
			pP[l] = P[l] + 1 + Sz2[0] + SzPl;

			pIm = Im + 1 + Sz2[0] + SzPl;
			OutIt.GoToBegin();

			for (int k = 0; k < Sz[2]; k++)
			{
				for (int j = 0; j < Sz[1]; j++)
				{
					for (int i = 0; i < Sz[0]; i++)
					{
						*pIm = (*pIm - Lam*(Sc[0]*(*pP[0] - *(pP[0] - 1)) + Sc[1]*(*pP[1] - *(pP[1] - Sz2[0])) + Sc[2]*(*pP[2] - *(pP[2] - SzPl))));
						for (int l = 0; l < Dim; l++)
						pP[l]++;
						OutIt.Set(static_cast<typename TOutputImage::PixelType>(*pIm++));
						++OutIt;
					}
					pIm += 2;
					for (int l = 0; l < Dim; l++)
					pP[l] += 2;
				}
				pIm += 2 * Sz2[0];
				for (int l = 0; l < Dim; l++)
				pP[l] += 2*Sz2[0];
			}
		}
		else
		{
		}
	}

	template<typename TInputImage, typename TOutputImage>
	void
	TotalVariationMinimization<TInputImage, TOutputImage>
	::EngineSliceBySlice()
	{
		uint32 Dim=TInputImage::ImageDimension;
		uint32 Sz2[3];
		uint32 SzPl;
		InPxType *pIm;
		typename TOutputImage::Pointer OutIm = this->GetOutput();
		ImageRegionConstIterator<TInputImage> InIt(const_cast< TInputImage* >(this->GetInput()),
				const_cast< TInputImage* >(this->GetInput())->GetRequestedRegion());
		ImageRegionIterator<TOutputImage> OutIt(this->GetOutput(), this->GetOutput()->GetRequestedRegion());

		if (Dim == 2)
		Sz[2] = 1;

		for (int i = 0; i<Dim; i++)
		Sz2[i] = Sz[i] + 2;

		OutIt.GoToBegin();
		InIt.GoToBegin();

		for (int k = 0; k < Sz[2]; k++)
		{
			pIm = Im + Sz2[0] + 1;
			for (int j = 0; j < Sz[1]; j++)
			{
				for (int i = 0; i < Sz[0]; i++)
				{
					*pIm++ = static_cast<InPxType>(InIt.Get());
					++InIt;
				}
				pIm += 2;
			}

			TVmin2D();

			pIm = Im + Sz2[0] + 1;
			for (int j = 0; j < Sz[1]; j++)
			{
				for (int i = 0; i < Sz[0]; i++)
				{
					OutIt.Set(static_cast<typename TOutputImage::PixelType>(*pIm++));
					++OutIt;
				}
				pIm += 2;
			}
		}
	}

	template<typename TInputImage, typename TOutputImage>
	void
	TotalVariationMinimization<TInputImage, TOutputImage>
	::PrintSelf(std::ostream & os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "Lambda: " << Lm << std::endl;
		os << indent << "Iteration Num: " << It << std::endl;
		os << indent << "To: " << It << std::endl;
		os << indent << "Iteration Num: " << It << std::endl;
	}

	template<typename TInputImage, typename TOutputImage>
	TotalVariationMinimization<TInputImage, TOutputImage>
	::~TotalVariationMinimization()
	{
		if (Allocated)
		{
			for (int i = 0; i<Dm; i++)
			delete[] P[i];
			delete[] Im;
			Allocated = false;
			TotalVox = 1;
		}
	}

}

#endif
