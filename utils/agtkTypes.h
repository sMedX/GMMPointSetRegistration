#ifndef __agtkTypes_h
#define __agtkTypes_h

#include <numeric>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkMesh.h>
#include <itkPoint.h>
#include <itkVectorContainer.h>
#include <itkImageMaskSpatialObject.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkPasteImageFilter.h>
#include <itkDefaultDynamicMeshTraits.h>
#include <itkIntTypes.h>
#include <itkShapeLabelObject.h>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>

namespace agtk
{
typedef signed char Int8Pixel;
typedef signed short Int16Pixel;
typedef signed long Int32Pixel;
typedef unsigned char UInt8Pixel;
typedef unsigned short UInt16Pixel;
typedef unsigned long UInt32Pixel;
typedef UInt8Pixel BinaryPixel;
typedef float FloatPixel;
typedef double DoublePixel;
typedef float RealPixel;

const unsigned int IMAGE_DIM_2 = 2;
const unsigned int IMAGE_DIM_3 = 3;

const BinaryPixel INSIDE_BINARY_VALUE = 1;
const BinaryPixel OUTSIDE_BINARY_VALUE = 0;

// ITK type aliases
typedef itk::Image<Int8Pixel, IMAGE_DIM_2> Int8Image2D;
typedef itk::Image<UInt8Pixel, IMAGE_DIM_2> UInt8Image2D;
typedef itk::Image<Int16Pixel, IMAGE_DIM_2> Int16Image2D;
typedef itk::Image<UInt16Pixel, IMAGE_DIM_2> UInt16Image2D;
typedef itk::Image<Int32Pixel, IMAGE_DIM_2> Int32Image2D;
typedef itk::Image<UInt32Pixel, IMAGE_DIM_2> UInt32Image2D;
typedef itk::Image<BinaryPixel, IMAGE_DIM_2> BinaryImage2D;
typedef itk::Image<FloatPixel, IMAGE_DIM_2> FloatImage2D;
typedef itk::Image<DoublePixel, IMAGE_DIM_2> DoubleImage2D;
typedef itk::Image<RealPixel, IMAGE_DIM_2> RealImage2D;

typedef itk::Image<Int8Pixel, IMAGE_DIM_3> Int8Image3D;
typedef itk::Image<UInt8Pixel, IMAGE_DIM_3> UInt8Image3D;
typedef itk::Image<Int16Pixel, IMAGE_DIM_3> Int16Image3D;
typedef itk::Image<UInt16Pixel, IMAGE_DIM_3> UInt16Image3D;
typedef itk::Image<Int32Pixel, IMAGE_DIM_3> Int32Image3D;
typedef itk::Image<UInt32Pixel, IMAGE_DIM_3> UInt32Image3D;
typedef itk::Image<BinaryPixel, IMAGE_DIM_3> BinaryImage3D;
typedef itk::Image<FloatPixel, IMAGE_DIM_3> FloatImage3D;
typedef itk::Image<DoublePixel, IMAGE_DIM_3> DoubleImage3D;
typedef itk::Image<RealPixel, IMAGE_DIM_3> RealImage3D;

typedef itk::VectorImage<Int8Pixel, IMAGE_DIM_2> Int8VectorImage2D;
typedef itk::VectorImage<UInt8Pixel, IMAGE_DIM_2> UInt8VectorImage2D;
typedef itk::VectorImage<Int16Pixel, IMAGE_DIM_2> Int16VectorImage2D;
typedef itk::VectorImage<UInt16Pixel, IMAGE_DIM_2> UInt16VectorImage2D;
typedef itk::VectorImage<Int32Pixel, IMAGE_DIM_2> Int32VectorImage2D;
typedef itk::VectorImage<UInt32Pixel, IMAGE_DIM_2> UInt32VectorImage2D;
typedef itk::VectorImage<BinaryPixel, IMAGE_DIM_2> BinaryVectorImage2D;
typedef itk::VectorImage<FloatPixel, IMAGE_DIM_2> FloatVectorImage2D;

typedef itk::VectorImage<Int8Pixel, IMAGE_DIM_3> Int8VectorImage3D;
typedef itk::VectorImage<UInt8Pixel, IMAGE_DIM_3> UInt8VectorImage3D;
typedef itk::VectorImage<Int16Pixel, IMAGE_DIM_3> Int16VectorImage3D;
typedef itk::VectorImage<UInt16Pixel, IMAGE_DIM_3> UInt16VectorImage3D;
typedef itk::VectorImage<Int32Pixel, IMAGE_DIM_3> Int32VectorImage3D;
typedef itk::VectorImage<UInt32Pixel, IMAGE_DIM_3> UInt32VectorImage3D;
typedef itk::VectorImage<BinaryPixel, IMAGE_DIM_3> BinaryVectorImage3D;
typedef itk::VectorImage<FloatPixel, IMAGE_DIM_3> FloatVectorImage3D;

typedef itk::FixedArray<float, IMAGE_DIM_3> FixedArrayPixelType3D;
typedef itk::Image<FixedArrayPixelType3D, IMAGE_DIM_3> FloatFixedArrayImage3D;

typedef itk::CovariantVector <float, IMAGE_DIM_3> FloatCovariantVector3D;
typedef itk::Image<FloatCovariantVector3D, IMAGE_DIM_3> FloatCovariantVectorImage3D;

typedef itk::ImageBase<IMAGE_DIM_2>::SizeType Image2DSize;
typedef itk::ImageBase<IMAGE_DIM_2>::IndexType Image2DIndex;
typedef itk::ImageBase<IMAGE_DIM_2>::OffsetType Image2DOffset;
typedef itk::ImageBase<IMAGE_DIM_2>::RegionType Image2DRegion;
typedef itk::ImageBase<IMAGE_DIM_2>::PointType Image2DPoint;
typedef itk::ImageBase<IMAGE_DIM_2>::SpacingType Image2DSpacing;

typedef itk::ImageBase<IMAGE_DIM_3>::SizeType Image3DSize;
typedef itk::ImageBase<IMAGE_DIM_3>::IndexType Image3DIndex;
typedef itk::ImageBase<IMAGE_DIM_3>::OffsetType Image3DOffset;
typedef itk::ImageBase<IMAGE_DIM_3>::RegionType Image3DRegion;
typedef itk::ImageBase<IMAGE_DIM_3>::PointType Image3DPoint;
typedef itk::ImageBase<IMAGE_DIM_3>::SpacingType Image3DSpacing;

typedef itk::Point<double, IMAGE_DIM_3> Point3D;

typedef itk::VectorContainer<unsigned int, Point3D> Point3DContainer;
typedef itk::VectorContainer<unsigned int, Image3DIndex> Image3DIndexContainer;

typedef itk::ImageMaskSpatialObject<IMAGE_DIM_3> ImageMaskSpatialObject3D;
typedef itk::BinaryBallStructuringElement<BinaryPixel, IMAGE_DIM_3> BinaryBallStructuringElement3D;

typedef itk::ShapeLabelObject<itk::SizeValueType, IMAGE_DIM_2> ShapeLabelObject2D;
typedef itk::ShapeLabelObject<itk::SizeValueType, IMAGE_DIM_3> ShapeLabelObject3D;

// ITK Filters
typedef itk::RegionOfInterestImageFilter<FloatImage3D, FloatImage3D> CropFloatImage3DFilter;
typedef itk::RegionOfInterestImageFilter<BinaryImage3D, BinaryImage3D> CropBinaryImage3DFilter;

typedef itk::PasteImageFilter<FloatImage3D> PasteFloatImage3DFilter;
typedef itk::PasteImageFilter<BinaryImage3D> PasteBinaryImage3DFilter;

// ITK meshes
typedef itk::DefaultDynamicMeshTraits<double, IMAGE_DIM_3, IMAGE_DIM_3, double, double> DoubleTriangleMeshTraits;
typedef itk::Mesh<double, IMAGE_DIM_3, DoubleTriangleMeshTraits> DoubleTriangleMesh3D;
typedef itk::Mesh<float, IMAGE_DIM_3> FloatTriangleMesh3D;

// VTK type aliases
typedef vtkSmartPointer<vtkImageData> VtkImageDataPointer;
typedef vtkSmartPointer<vtkPolyData> VtkPolyDataPointer;
typedef vtkSmartPointer<vtkIdList> VtkIdListPointer;

// numeric limits
typedef std::numeric_limits<double> DoubleLimits;
typedef std::numeric_limits<float> FloatLimits;
typedef std::numeric_limits<BinaryPixel> BinaryLimits;
}

#endif // __agtkTypes_h
