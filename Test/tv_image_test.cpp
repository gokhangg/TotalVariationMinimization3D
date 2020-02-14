
#include "tv_image.h"
#include "gtest/gtest.h"
#include <vector>
#include <cmath>


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

constexpr float EPSILON = 1e-2;

auto GetInitializedData(const float diff = 1.f)
{
	std::vector<unsigned int> imSize{ 11, 12, 13 };
	float in = 0.f;
	auto fillVal = [&](const float&)
	{
		in += diff;
		return in;
	};
	TVimage<> tvImage(imSize);
	tvImage.transform(fillVal);
	return tvImage;
}

TEST(TVimage, Gradient)
{
	const float diff = 10.99f;
	auto tvImage = GetInitializedData(diff);
	auto grad = tvImage.getGradient();

	auto stride = tvImage.getStride();

	//Due to precision error, EPSILON should be applied;
	EXPECT_NEAR(grad[0][0], diff, EPSILON);
	EXPECT_NEAR(grad[1][0], stride[1] * diff, EPSILON);
	EXPECT_NEAR(grad[2][10], stride[2] * diff, EPSILON);
}


TEST(TVimage, Divergence)
{
	const float diff = 100.99f;
	auto tvImage = GetInitializedData(diff);
	auto grad = tvImage.getGradient();
	auto div = TVimage<>::getDivergence(grad);

	auto stride = tvImage.getStride();

	//Due to precision error, EPSILON should be applied;
	EXPECT_NEAR(div[0], 0, EPSILON);
	EXPECT_NEAR(div[stride[1]], 0, EPSILON);
	EXPECT_NEAR(div[stride[2]], 0, EPSILON);
}


//Test for operators should be added.