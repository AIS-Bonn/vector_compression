// Unit test for vector_compression

#include <vector_compression/vector_compression.h>

#include "../src/contrib/angles.h"

#include <iostream>

#include "catch.hpp"


TEST_CASE("weird quat", "[vector_compression]")
{
	Eigen::Quaternionf source(-0.414, 0.0, 0.0, 0.910);

	uint8_t buffer[5];

	vector_compression::encodeQuaternion(source, buffer);

	Eigen::Quaternionf out = vector_compression::decodeQuaternion(buffer);

	INFO("In: " << source.coeffs().transpose() << ", Out: " << out.coeffs().transpose());

	Eigen::Quaternionf diff = out * source.inverse();

	Eigen::AngleAxisf aa;
	aa = diff;

	float angleDiff = angles::normalize_angle(aa.angle());
	REQUIRE(fabsf(angleDiff) < 0.5f * M_PI / 180.0f);
}

TEST_CASE("Identity", "[vector_compression]")
{
	Eigen::Quaternionf quat = Eigen::Quaternionf::Identity();

	uint8_t buffer[5];
	vector_compression::encodeQuaternion(quat, buffer);

	Eigen::Quaternionf out = vector_compression::decodeQuaternion(buffer);

	Eigen::AngleAxisf aa;
	aa = out;

	REQUIRE(aa.angle() == 0);
}

TEST_CASE("Simple axis", "[vector_compression]")
{
	Eigen::Vector3f axis(1.0, 1.0, 1.0);
	axis.normalize();

	float max_err = 0.0;

	for(int i = -100; i < 100; ++i)
	{
		uint8_t buffer[5];

		Eigen::Quaternionf quat;
		quat = Eigen::AngleAxisf(M_PI / 100 * i, axis);

		INFO("Input: " << (M_PI / 100 * i));

		vector_compression::encodeQuaternion(quat, buffer);

		Eigen::Quaternionf out = vector_compression::decodeQuaternion(buffer);

		Eigen::AngleAxisf outa;
		outa = out;
		INFO("Output: " << outa.angle());

		Eigen::Quaternionf diff = out * quat.inverse();
		Eigen::AngleAxisf aa;
		aa = diff;

		float angleDiff = angles::normalize_angle(aa.angle());
		REQUIRE(fabsf(angleDiff) < 0.5);

		max_err = std::max(max_err, fabsf(angleDiff));
	}

	printf("max error: %f (%f deg)\n", max_err, max_err * 180.0 / M_PI);
}

TEST_CASE("Z axis", "[vector_compression]")
{
	Eigen::Vector3f axis(0.0, 0.0, 1.0);
	axis.normalize();

	float max_err = 0.0;

	for(int i = -100; i < 100; ++i)
	{
		uint8_t buffer[5];

		Eigen::Quaternionf quat;
		quat = Eigen::AngleAxisf(M_PI / 100 * i, axis);

		INFO("Input: " << (M_PI / 100 * i));

		vector_compression::encodeQuaternion(quat, buffer);

		Eigen::Quaternionf out = vector_compression::decodeQuaternion(buffer);

		Eigen::AngleAxisf outa;
		outa = out;
		INFO("Output: " << outa.angle());

		Eigen::Quaternionf diff = out * quat.inverse();
		Eigen::AngleAxisf aa;
		aa = diff;

		float angleDiff = angles::normalize_angle(aa.angle());
		REQUIRE(fabsf(angleDiff) < 0.5);

		max_err = std::max(max_err, fabsf(angleDiff));
	}

	printf("Z axis: max error: %f (%f deg)\n", max_err, max_err * 180.0 / M_PI);
}

TEST_CASE("Lattice", "[vector_compression]")
{
	float max_err = 0;

	for(int i = -100; i < -50; ++i)
	{
		for(int j = -100; j < -50; ++j)
		{
			for(int k = -100; k < -50; ++k)
			{
				Eigen::Vector3f vec(i * 0.01, j * 0.01, k * 0.01);

				Eigen::Vector3i lattice = vector_compression::vectorToLattice(vec, 0.1);

				Eigen::Vector3f backProj = vector_compression::latticeToVector(lattice, 0.1);

				float dist = (vec - backProj).norm();

				INFO("vector: " << vec.transpose());
				INFO("lattice: " << lattice.transpose());
				INFO("backProj: " << backProj.transpose());
				REQUIRE(dist < 0.1 + 1e-5);

				max_err = std::max(max_err, dist);
			}
		}
	}

	printf("lattice: max_err = %f\n", max_err);
}

TEST_CASE("Aligned lattice", "[vector_compression]")
{
	float max_err = 0;
	float max_simpleErr = 0;
	Eigen::Vector3f maxCoeff = Eigen::Vector3f::Zero();

	const float r = 1.0f / pow(2*1000, 1.0f/3.0f);

	for(int i = -100; i < -50; ++i)
	{
		for(int j = -100; j < -50; ++j)
		{
			for(int k = -100; k < -50; ++k)
			{
				Eigen::Vector3f vec(i * 0.01, j * 0.01, k * 0.01);

				Eigen::Vector3i lattice = vector_compression::vectorToAlignedLattice(vec, r);

				Eigen::Vector3f backProj = vector_compression::alignedLatticeToVector(lattice, r);

				float dist = (vec - backProj).norm();

				INFO("vector: " << vec.transpose());
				INFO("lattice: " << lattice.transpose());
				INFO("backProj: " << backProj.transpose());
				REQUIRE(dist < sqrtf(2.0)*r + 1e-5);

				Eigen::Vector3i simpleLattice;
				simpleLattice <<
					std::round(vec.x() * 10),
					std::round(vec.y() * 10),
					std::round(vec.z() * 10)
				;

				Eigen::Vector3f simpleBackProj = simpleLattice.cast<float>() / 10.0f;
				float simpleDist = (vec - simpleBackProj).norm();

				max_err = std::max(max_err, dist);
				max_simpleErr = std::max(max_simpleErr, simpleDist);
				maxCoeff.x() = std::max(maxCoeff.x(), fabsf(lattice.x()));
				maxCoeff.y() = std::max(maxCoeff.y(), fabsf(lattice.y()));
				maxCoeff.z() = std::max(maxCoeff.z(), fabsf(lattice.z()));
			}
		}
	}

	printf("lattice: max_err = %f, maxCoeff: %f, %f, %f (in comparison to simple err %f)\n", max_err, maxCoeff.x(), maxCoeff.y(), maxCoeff.z(), max_simpleErr);
}

TEST_CASE("Vector", "[vector_compression]")
{
	float max_err = 0;

	for(int i = -100; i < -50; ++i)
	{
		for(int j = -100; j < -50; ++j)
		{
			for(int k = -100; k < -50; ++k)
			{
				uint8_t buf[6];

				Eigen::Vector3f vec(i * 0.01, j * 0.01, k * 0.01);

				vector_compression::encode3DVector6(vec, 2.0, buf);

				Eigen::Vector3f backProj = vector_compression::decode3DVector6(2.0, buf);

				float dist = (vec - backProj).norm();

				max_err = std::max(max_err, dist);
			}
		}
	}

	printf("vector compression: max_err = %f\n", max_err);
}

TEST_CASE("VectorCompression", "[vector_compression]")
{
	float max_err = 0;
	vector_compression::VectorCompression compression(24);

	for(int i = -100; i < 100; i += 4)
	{
		for(int j = -100; j < -50; ++j)
		{
			for(int k = -100; k < -50; ++k)
			{
				uint64_t buf;

				Eigen::Vector3f vec(i * 0.01, j * 0.007, k * 0.004);

				buf = compression.compress(vec);

				REQUIRE(buf < (1 << 24));

				Eigen::Vector3f backProj = compression.uncompress(buf);

				float dist = (vec - backProj).norm();

				REQUIRE(dist < compression.guaranteedPrecision());

				max_err = std::max(max_err, dist);
			}
		}
	}

	// Test that it's better than simple cubic packing
	float simpleCubicPrec = sqrt(3.0) * pow((2.0*2.0*2.0) / (1LL << 24), 1.0f/3.0f);
	REQUIRE(max_err < simpleCubicPrec);

	printf("vector compression: max_err = %f, guaranteed = %f (in comp. to %f for cubic packing)\n", max_err, compression.guaranteedPrecision(), simpleCubicPrec);
}

TEST_CASE("VectorCompression with rectangular shape", "[vector_compression]")
{
	float max_err = 0;
	vector_compression::VectorCompression compression(24);

	compression.setExtent(Eigen::Vector3f(3.0, 1.0, 1.0));

	for(int i = -100; i < 100; i += 4)
	{
		for(int j = -100; j < -50; ++j)
		{
			for(int k = -100; k < -50; ++k)
			{
				uint64_t buf;

				Eigen::Vector3f vec(3.0 * i * 0.01, j * 0.007, k * 0.004);

				buf = compression.compress(vec);

				REQUIRE(buf < (1 << 24));

				Eigen::Vector3f backProj = compression.uncompress(buf);

				float dist = (vec - backProj).norm();

				max_err = std::max(max_err, dist);
			}
		}
	}

	// Test that it's better than simple cubic packing
	float simpleCubicPrec = sqrt(3.0) * pow((6.0*2.0*2.0) / (1LL << 24), 1.0f/3.0f);
	REQUIRE(max_err < simpleCubicPrec);

	printf("vector compression: max_err = %f (in comp. to %f)\n", max_err, simpleCubicPrec);
}
