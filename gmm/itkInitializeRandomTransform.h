#pragma once

#include <itkObject.h>
#include <itkTranslationTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkScaleSkewVersor3DTransform.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

namespace itk
{
  template<typename TParametersValueType>
  class InitializeRandomTransform : public Object
  {
  public:
    /** Standard class typedefs. */
    typedef InitializeRandomTransform  Self;
    typedef Object                     Superclass;
    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    itkTypeMacro(InitializeRandomTransform, Object);

    enum Transform
    {
      Translation = 0,
      Versor3D = 1,
      Similarity = 2,
      ScaleSkewVersor3D = 3
    };

    /** typedefs */
    itkStaticConstMacro(PointDimension, unsigned int, 3U);
    static_assert(PointDimension == 3U, "Invalid dimension. Dimension 3 is supported.");

    typedef TParametersValueType                                           ParametersValueType;
    typedef typename itk::Transform<ParametersValueType, PointDimension>   TransformType;
    typedef typename TransformType::InputPointType                         InputPointType;
    typedef typename TransformType::OutputPointType                        OutputPointType;
    typedef typename TransformType::OutputVectorType                       OutputVectorType;

    typedef typename TransformType::ParametersType                         ParametersType;
    typedef Versor<ParametersValueType>                                    VersorType;
    typedef Vector<ParametersValueType, PointDimension>                    AxisType;

    /** Type of the transform generator. */
    typedef itk::Statistics::MersenneTwisterRandomVariateGenerator TransformGeneratorType;

    // Get transform
    itkGetObjectMacro(Transform, TransformType);

    // Set/Get type of transform
    itkSetEnumMacro(TypeOfTransform, Transform);
    itkGetEnumMacro(TypeOfTransform, Transform);
    void SetTypeOfTransform(const size_t & type) { this->SetTypeOfTransform(static_cast<Transform>(type)); }

    // Set/Get the center 
    itkSetMacro(Center, InputPointType);
    itkGetMacro(Center, InputPointType);

    // Set/Get the bounds
    itkSetMacro(TranslationBounds, double);
    itkGetMacro(TranslationBounds, double);

    itkSetMacro(RotationBounds, double);
    itkGetMacro(RotationBounds, double);

    itkSetMacro(ScalingBounds, double);
    itkGetMacro(ScalingBounds, double);

    itkSetMacro(SkewBounds, double);
    itkGetMacro(SkewBounds, double);

    void Initialize()
    {
      switch (m_TypeOfTransform) {
      case Transform::Translation: {
        // Translation transform
        typedef typename itk::TranslationTransform<TParametersValueType, PointDimension> TranslationTransformType;
        typename TranslationTransformType::Pointer transform = TranslationTransformType::New();

        m_NumberOfTranslationComponents = 3;

        ParametersType parameters = transform->GetParameters();

        for (size_t n = 0; n < m_NumberOfTranslationComponents; ++n)
        {
          parameters[n] = m_TransformGenerator->GetUniformVariate(-m_TranslationBounds, m_TranslationBounds);
        }

        m_Transform = transform;
        m_Transform->SetParameters(parameters);
        break;
      }
      case Transform::Versor3D: {
        // VersorRigid3DTransform
        typedef itk::VersorRigid3DTransform<TParametersValueType> VersorRigid3DTransformType;
        typename VersorRigid3DTransformType::Pointer transform = VersorRigid3DTransformType::New();
        transform->SetIdentity();
        transform->SetCenter(m_Center);

        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;

        ParametersType parameters = transform->GetParameters();

        int count = 0;

        // rotation
        this->InitializeVersor();

        for (size_t n = 0; n < m_NumberOfRotationComponents; ++n, ++count) 
        {
          parameters[count] = m_Versor.GetVnlQuaternion()[n];
        }

        // translation
        for (size_t n = 0; n < m_NumberOfTranslationComponents; ++n, ++count)
        {
          parameters[count] = m_TransformGenerator->GetUniformVariate(-m_TranslationBounds, m_TranslationBounds);
        }

        m_Transform = transform;
        m_Transform->SetParameters(parameters);
        break;
      }

      case Transform::Similarity:{
        // Similarity3DTransform
        typedef itk::Similarity3DTransform<TParametersValueType> Similarity3DTransformType;
        typename Similarity3DTransformType::Pointer transform = Similarity3DTransformType::New();
        transform->SetIdentity();
        transform->SetCenter(m_Center);

        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;
        m_NumberOfScalingComponents = 1;

        ParametersType parameters = transform->GetParameters();

        int count = 0;

        // rotation
        this->InitializeVersor();

        for (size_t n = 0; n < m_NumberOfRotationComponents; ++n, ++count) 
        {
          parameters[count] = m_Versor.GetVnlQuaternion()[n];
        }

        // translation
        for (size_t n = 0; n < m_NumberOfTranslationComponents; ++n, ++count) 
        {
          parameters[count] = m_TransformGenerator->GetUniformVariate(-m_TranslationBounds, m_TranslationBounds);
        }

        // scaling
        for (size_t n = 0; n < m_NumberOfScalingComponents; ++n, ++count)
        {
          parameters[count] = m_TransformGenerator->GetUniformVariate(1 - m_ScalingBounds, 1 + m_ScalingBounds);
        }

        m_Transform = transform;
        m_Transform->SetParameters(parameters);
        break;
      }

      case Transform::ScaleSkewVersor3D:{
        typedef itk::ScaleSkewVersor3DTransform<TParametersValueType> ScaleSkewVersor3DTransformType;
        typename ScaleSkewVersor3DTransformType::Pointer transform = ScaleSkewVersor3DTransformType::New();
        transform->SetIdentity();
        transform->SetCenter(m_Center);

        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;
        m_NumberOfScalingComponents = 3;
        m_NumberOfSkewComponents = 6;

        m_Transform = transform;
        ParametersType parameters = m_Transform->GetParameters();

        int count = 0;

        // rotation
        this->InitializeVersor();

        for (size_t n = 0; n < m_NumberOfRotationComponents; ++n, ++count) 
        {
          parameters[count] = m_Versor.GetVnlQuaternion()[n];
        }

        // translation
        for (size_t n = 0; n < m_NumberOfTranslationComponents; ++n, ++count) 
        {
          parameters[count] = m_TransformGenerator->GetUniformVariate(-m_TranslationBounds, m_TranslationBounds);
        }

        // scaling
        for (size_t n = 0; n < m_NumberOfScalingComponents; ++n, ++count) 
        {
          parameters[count] = m_TransformGenerator->GetUniformVariate(1 - m_ScalingBounds, 1 + m_ScalingBounds);
        }

        // skew
        for (size_t n = 0; n < m_NumberOfSkewComponents; ++n, ++count) 
        {
          parameters[count] = m_TransformGenerator->GetUniformVariate(-m_SkewBounds, m_SkewBounds);
        }

        m_Transform->SetParameters(parameters);
        break;
      }
      default:
        itkExceptionMacro(<< "Invalid type of the input transform");
      }
    }

  protected:
    InitializeRandomTransform() 
    {
      m_TransformGenerator = TransformGeneratorType::New();
      m_TransformGenerator->Initialize(m_RandomSeed);
    }
    virtual ~InitializeRandomTransform() {};

    void InitializeVersor() 
    {
      m_Angle = m_TransformGenerator->GetUniformVariate(-m_RotationBounds, m_RotationBounds);

      for (size_t n = 0; n < PointDimension; ++n) {
        m_Axis[n] = m_TransformGenerator->GetUniformVariate(-1, 1);
      }
      m_Versor.Set(m_Axis, m_Angle);
    }

    virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE
    {
      Superclass::PrintSelf(os, indent);
      os << std::endl;

      os << "Transform " << m_Transform->GetTransformTypeAsString() << std::endl;
      os << std::endl;

      os << "Bounds" << std::endl;
      os << indent << "Translation bounds " << m_TranslationBounds << std::endl;
      os << indent << "Rotation bounds    " << m_RotationBounds << std::endl;
      os << indent << "Scaling bounds     " << m_ScalingBounds << std::endl;
      os << indent << "Skew bounds        " << m_SkewBounds << std::endl;
      os << std::endl;

      os << "Parameters" << std::endl;
      os << indent << "Center           " << m_Center << std::endl;
      os << indent << "Fixed parameters " << m_Transform->GetFixedParameters() << ", " << m_Transform->GetNumberOfFixedParameters() << std::endl;
      os << indent << "Parameters       " << m_Transform->GetParameters() << ", " << m_Transform->GetNumberOfParameters() << std::endl;
      os << indent << "Versor           " << m_Versor << ", angle " << m_Versor.GetAngle() <<  ", axis " << m_Versor.GetAxis() << std::endl;
      os << std::endl;
    }

    Transform m_TypeOfTransform = Transform::Similarity;
    typename TransformType::Pointer m_Transform = nullptr;

    InputPointType m_Center;
    VersorType m_Versor;
    AxisType m_Axis;
    double m_Angle;

    typename TransformGeneratorType::Pointer m_TransformGenerator;
    int m_RandomSeed = 0;

    double m_TranslationBounds = 100;
    double m_RotationBounds = itk::Math::pi_over_4;
    double m_ScalingBounds = 0.1;
    double m_SkewBounds = 0.1;

    size_t m_NumberOfTranslationComponents = 0;
    size_t m_NumberOfRotationComponents = 0;
    size_t m_NumberOfScalingComponents = 0;
    size_t m_NumberOfSkewComponents = 0;

  private:
    InitializeRandomTransform(const Self &) ITK_DELETE_FUNCTION;
    void operator=(const Self &) ITK_DELETE_FUNCTION;
  };
}
