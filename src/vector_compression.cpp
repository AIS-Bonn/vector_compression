// Compression methods for vectors
// Author: Max Schwarz <max.schwarz@uni-bonn.de>

#include <vector_compression/vector_compression.h>

#include <Eigen/Core>

#include "contrib/angles.h"

#include <iostream>

#define DEBUG 0

namespace vector_compression
{

// The "hemi-oct" compression scheme is taken from
//  Cigolle, Zina H., Sam Donow, and Daniel Evangelakos.
//  "A survey of efficient representations for independent unit vectors."
//  Journal of Computer Graphics Techniques Vol 3.2 (2014).

// Assume normalized input on +Z hemisphere.
// Output is [-1, 1]
static Eigen::Vector2f float32x3_to_hemioct(const Eigen::Vector3f& v)
{
	// Project the hemisphere onto the hemi-octahedron,
	// and then project into the xy plane
	Eigen::Vector2f p = v.head<2>() * (1.0f / (fabsf(v.x()) + fabsf(v.y()) + v.z()));

	// Rotate and scale the center diamond to the unit square
	return Eigen::Vector2f(p.x() + p.y(), p.x() - p.y());
}

static Eigen::Vector3f hemioct_to_float32x3(const Eigen::Vector2f& e)
{
	// Rotate and scale the unit square back to the center diamond
	Eigen::Vector2f temp = 0.5f * Eigen::Vector2f(e.x() + e.y(), e.x() - e.y());
	Eigen::Vector3f v;
	v << temp, 1.0 - fabsf(temp.x()) - fabsf(temp.y());

	return v.normalized();
}

int32_t float_to_snormb(float v, int b)
{
	return std::round(std::min(1.0f, std::max(-1.0f, v)) * ((1 << (b-1)) - 1));
}

float snormb_to_float(int32_t v, int b)
{
	// Sign extend
	if(v & (1 << (b-1)))
		v |= 0xFFFFFFFF & ~((1 << b) - 1);

	return std::min(1.0f, std::max(-1.0f, ((float)v) / ((1 << (b-1)) - 1)));
}

void encodeQuaternion(const Eigen::Quaternionf& quat, uint8_t* dest)
{
	// Recover rotation angle & axis
	float angle = 2.0f * acosf(quat.w());
	angle = angles::normalize_angle(angle);

	float norm = sqrtf(1.0f - quat.w()*quat.w());

	Eigen::Vector3f axis;

	if(norm != 0)
	{
		axis << quat.x() / norm, quat.y() / norm, quat.z() / norm;
	}
	else
		axis << 0.0f, 0.0f, 0.0f;

	// Map to the upper hemisphere
	if(axis.z() < 0)
	{
		axis *= -1.0f;
		angle = -angle;
	}

	Eigen::Vector2f hemi = float32x3_to_hemioct(axis);

	int32_t hemi_a = float_to_snormb(hemi.x(), 12);
	int32_t hemi_b = float_to_snormb(hemi.y(), 12);

	int32_t angle_sn = float_to_snormb(angle / M_PI, 16);

#if DEBUG
	printf("encode: %f, %f, %f, angle_sn: %u (0x%X)\n", hemi.x(), hemi.y(), angle, angle_sn, angle_sn);
	printf("angle scale: %d\n", ((1 << (16-1)) - 1));
#endif

	dest[0] = hemi_a;
	dest[1] = ((hemi_a >> 8) & 0xF) | ((hemi_b & 0xF) << 4);
	dest[2] = hemi_b >> 4;
	dest[3] = angle_sn;
	dest[4] = angle_sn >> 8;
}

Eigen::Quaternionf decodeQuaternion(const uint8_t* src)
{
	int32_t hemi_a = src[0] | (((uint32_t)src[1] & 0xF) << 8);
	int32_t hemi_b = (src[1] >> 4) | (((uint32_t)src[2]) << 4);
	int32_t angle_sn = src[3] | (((uint32_t)src[4]) << 8);

	Eigen::Vector2f hemi(
		snormb_to_float(hemi_a, 12), snormb_to_float(hemi_b, 12)
	);

	Eigen::Vector3f axis = hemioct_to_float32x3(hemi);

	float angle = snormb_to_float(angle_sn, 16) * M_PI;

#if DEBUG
	printf("decode: %f, %f, %f, angle_sn: %d\n", hemi.x(), hemi.y(), angle, angle_sn);
#endif

	return Eigen::Quaternionf(Eigen::AngleAxisf(angle, axis));
}

Eigen::Vector3f latticeToVector(const Eigen::Vector3i& lattice, float r)
{
	return r * Eigen::Vector3f(
		lattice.y() + lattice.z(),
		lattice.x() + lattice.z(),
		lattice.x() + lattice.y()
	);
}

Eigen::Vector3f alignedLatticeToVector(const Eigen::Vector3i& lattice, float r)
{
	Eigen::Vector3f m = lattice.cast<float>();

	Eigen::Vector3f x = m + m.y() * Eigen::Vector3f::UnitY() + ((lattice.x() + lattice.z()) % 2) * Eigen::Vector3f::UnitY();

	x *= r;

	return x;
}

static inline int sgn(float x)
{
	if(x >= 0)
		return 1;
	else
		return -1;
}

static inline float roundFloat(float x)
{
	return std::round(x);
}

Eigen::Vector3i vectorToAlignedLattice(const Eigen::Vector3f& p, float r)
{
	Eigen::Vector3f m = p / r;

	Eigen::Vector3f cartesian_lattice = m.unaryExpr(std::ptr_fun(roundFloat));
	Eigen::Vector3i lattice_int = cartesian_lattice.cast<int>();

	if(lattice_int.sum() % 2 != 0)
	{
		Eigen::Vector3f local = m - cartesian_lattice;
		int c;
		local.array().abs().maxCoeff(&c);

		lattice_int[c] += sgn(local[c]);
	}

	lattice_int.y() = (lattice_int.y() - ((lattice_int.x() + lattice_int.z()) % 2)) / 2;

	return lattice_int;
}

Eigen::Vector3i vectorToLattice(const Eigen::Vector3f& p, float r)
{
	Eigen::Vector3f normalized = p / r;

	Eigen::Vector3f cartesian_lattice = normalized.unaryExpr(std::ptr_fun(roundFloat));
	Eigen::Vector3i lattice_int = cartesian_lattice.cast<int>();

	if(lattice_int.sum() % 2 != 0)
	{
		Eigen::Vector3f local = normalized - cartesian_lattice;
		int c;
		local.array().abs().maxCoeff(&c);

		lattice_int[c] += sgn(local[c]);
	}

	Eigen::Matrix3i inv;
	inv <<
		-1, 1, 1,
		1, -1, 1,
		1, 1, -1
	;

	return (inv * lattice_int) / 2;
}

void encode3DVector6(const Eigen::Vector3f& vec, float extent, uint8_t* dest)
{
	// We are using 16 bit resolution per lattice index.
	float lattice_r = extent / (2.0f * ((1 << 15)-1));

	Eigen::Vector3i lattice = vectorToLattice(vec, lattice_r);

	for(int i = 0; i < 3; ++i)
	{
		lattice[i] = std::min(32767, std::max(-32768, lattice[i]));
		dest[2*i] = lattice[i];
		dest[2*i+1] = lattice[i] >> 8;
	}
}

Eigen::Vector3f decode3DVector6(float extent, const uint8_t* src)
{
	// We are using 16 bit resolution per lattice index.
	float lattice_r = extent / (2.0f * ((1 << 15)-1));

	Eigen::Vector3i lattice;
	for(int i = 0; i < 3; ++i)
	{
		int16_t idx = src[2*i] | (((uint16_t)src[2*i+1]) << 8);
		lattice[i] = idx;
	}

	return latticeToVector(lattice, lattice_r);
}

VectorCompression::VectorCompression(unsigned int bits)
 : m_bits(bits)
{
	setExtent(Eigen::Vector3f(1.0, 1.0, 1.0));
}

void VectorCompression::setExtent(const Eigen::Vector3f& extent)
{
	m_extent = extent;

	m_r = powf(sqrtf(2.0) * extent.x() * extent.y() * extent.z() / (1LL << m_bits), 1.0f / 3.0f) * sqrtf(2.0);

	m_latticeExtent.x() = std::floor(extent.x() / m_r);
	m_latticeExtent.y() = std::floor(extent.y() / m_r / 2.0f);
	m_latticeExtent.z() = std::floor(extent.z() / m_r);

	uint64_t high = 8 * m_latticeExtent.x() * m_latticeExtent.y() * m_latticeExtent.z();
	if(high > (1ULL << m_bits))
	{
		std::cerr << "lattice overflow: " << high << " > " << (1LL << 24) << "\n";
		std::cerr << "calculated radius: " << m_r;
		throw std::runtime_error("lattice overflow");
	}
}

uint64_t VectorCompression::compress(const Eigen::Vector3f& vec) const
{
	Eigen::Vector3i lattice = vectorToAlignedLattice(vec, m_r);

	// Shift to [0, 2*m_latticeExtent]
	lattice += m_latticeExtent;

	// Constrain
	lattice.x() = std::max(0, std::min(2*m_latticeExtent.x(), lattice.x()));
	lattice.y() = std::max(0, std::min(2*m_latticeExtent.y(), lattice.y()));
	lattice.z() = std::max(0, std::min(2*m_latticeExtent.z(), lattice.z()));

	uint64_t index = lattice.x() + lattice.y() * 2 * m_latticeExtent.x() + lattice.z() * 4 * m_latticeExtent.x() * m_latticeExtent.y();

	// Sanity check
	if(index >= (1ULL << m_bits))
	{
		std::cerr << "lattice overflow while compressing " << vec.transpose() << "\n";
		std::cerr << "radius: " << m_r << "\n";
		std::cerr << "lattice extents: " << m_latticeExtent.transpose() << "\n";
		std::cerr << "lattice coords: " << lattice.transpose() << "\n";
		std::cerr << "=> index: " << index << "\n";
		index = (1 << m_bits) - 1;
	}

	return index;
}

Eigen::Vector3f VectorCompression::uncompress(uint64_t compressed) const
{
	Eigen::Vector3i lattice;

	lattice.z() = compressed / (4 * m_latticeExtent.x() * m_latticeExtent.y());

	uint64_t remainder = compressed % (4 * m_latticeExtent.x() * m_latticeExtent.y());

	lattice.y() = remainder / (2*m_latticeExtent.x());
	lattice.x() = remainder % (2*m_latticeExtent.x());

	// Shift to [-extent, extent]
	lattice -= m_latticeExtent;

	return alignedLatticeToVector(lattice, m_r);
}

float VectorCompression::guaranteedPrecision() const
{
	return sqrt(2.0f) * m_r;
}

}
