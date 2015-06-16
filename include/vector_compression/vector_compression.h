// Compression methods for vectors
// Author: Max Schwarz <max.schwarz@uni-bonn.de>

#ifndef VECTOR_COMPRESSION_H
#define VECTOR_COMPRESSION_H

#include <Eigen/Geometry>
#include <Eigen/Core>

#include <stdint.h>

namespace vector_compression
{

//! @name snormb
//@{

/**
 * @brief Normalized signed float to binary
 *
 * @param v Input float in range [-1.0, 1.0]
 * @param b Output width in bits (including sign bit)
 * @return Binary representation of @a v in @a b bits
 **/
int32_t float_to_snormb(float v, int b);

/**
 * @brief Binary to normalized signed float
 *
 * @param v Binary representation in @a b bits
 * @param b Input width in bits (including sign bit)
 * @return Float in range [-1.0, 1.0]
 **/
float snormb_to_float(int32_t v, int b);
//@}

//! @name Quaternion compression
//@{

/**
 * @brief Quaternion to 5 bytes
 *
 * @param quat Input quaternion
 * @param dest Output array of size 5
 */
void encodeQuaternion(const Eigen::Quaternionf& quat, uint8_t* dest);

/**
 * @brief 5 Bytes to Quaternion
 *
 * @param src Input array of size 5
 * @return Quaternion
 **/
Eigen::Quaternionf decodeQuaternion(const uint8_t* src);
//@}

//! @name 3D vector compression
//@{

/**
 * @brief Rotated lattice to 3D vector
 *
 * This function maps lattice coordinates to a cubic space, which has its
 * vertices on the coordinate system axes.
 **/
Eigen::Vector3f latticeToVector(const Eigen::Vector3i& lattice, float r);
Eigen::Vector3i vectorToLattice(const Eigen::Vector3f& p, float r);

Eigen::Vector3f alignedLatticeToVector(const Eigen::Vector3i& lattice, float r);
Eigen::Vector3i vectorToAlignedLattice(const Eigen::Vector3f& p, float r);

void encode3DVector6(const Eigen::Vector3f& vec, float extent, uint8_t* dest);
Eigen::Vector3f decode3DVector6(float extent, const uint8_t* src);

/**
 * Flexible 3D vector compression / quantization using an FCC
 * (face-centered cubic) lattice. This approach has lower error than the
 * naive cubic packing (downsampling each dimension separately), since in that
 * case the error is much higher in diagonal directions.
 **/
class VectorCompression
{
public:
	/**
	 * Constructor.
	 *
	 * @param bits The number of bits to use for compression. Should be lower
	 *    than 64.
	 **/
	VectorCompression(unsigned int bits);

	//! @brief Number of bits per 3D vector
	inline unsigned int bits() const
	{ return m_bits; }

	/**
	 * @brief Set the covered extent.
	 *
	 * VectorCompression will compress points from -extent to extent with
	 * uniform precision. Points from outside the extent will be clamped
	 * to the extent range.
	 **/
	void setExtent(const Eigen::Vector3f& extent);

	/**
	 * @brief Compression
	 *
	 * @param point Input point in [-extent, extent]
	 * @return Integer with the lower bits() bits populated
	 **/
	uint64_t compress(const Eigen::Vector3f& point) const;

	/**
	 * @brief Uncompression
	 *
	 * @param compressed Integer with the lower bits() bits populated
	 * @return Uncompressed vector
	 **/
	Eigen::Vector3f uncompress(uint64_t compressed) const;

	/**
	 * @brief Precision
	 *
	 * @return The guaranteed precision (i.e. maximum error in arbitrary
	 *         direction)
	 **/
	float guaranteedPrecision() const;
private:
	unsigned int m_bits;
	Eigen::Vector3f m_extent;
	float m_r;

	Eigen::Vector3i m_latticeExtent;
};

//@}

}

#endif
