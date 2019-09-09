/* class of neural network where activation 
 function is ReLU on hidden and tanh on output */

class neuralnetwork
{
  ArrayList<double[][]> Weights = new ArrayList<double[][]>();
  ArrayList<double[]> Biases = new ArrayList<double[]>();
  ArrayList<double[]> Layers = new ArrayList<double[]>();
  ArrayList<Integer> LayersDimensions;

  neuralnetwork(ArrayList<Integer> layersDimensions)
  {
    NeuralNetworkFromArrayList(layersDimensions);
  }

  neuralnetwork(int input, int hiddenLayers, int hiddenNodes, int output)
  {
    ArrayList<Integer> layersDimensions = new ArrayList<Integer>();
    layersDimensions.add(input);
    for (int i=0; i<hiddenLayers; i++) 
    {
      layersDimensions.add(hiddenNodes);
    }
    layersDimensions.add(output);
    NeuralNetworkFromArrayList(layersDimensions);
  }

  private void NeuralNetworkFromArrayList(ArrayList<Integer> layersDimensions)
  {
    LayersDimensions = layersDimensions;
    Layers.add(new double[layersDimensions.get(0)]);
    for (int i=0; i<LayersDimensions.size()-1; i++)
    {
      Weights.add(new double[layersDimensions.get(i)][layersDimensions.get(i+1)]);
      Biases.add(new double[layersDimensions.get(i+1)]);
      Layers.add(new double[layersDimensions.get(i)]);
    }
  }
  

  double[] ComputeOutput(double[] input)
  {
    if (Layers.get(0).length == input.length)
    {
      Layers.set(0, input);
      for (int l = 1; l< Layers.size(); l++)
      {
        Layers.set(l, Mult(Weights.get(l-1), Layers.get(l-1))); //multiply by weights
        Sum(Layers.get(l), Biases.get(l-1)); //add bias
        if (l != Layers.size()-1) { //last layer
          ReLU(Layers.get(l)); //apply activation function
        } else
        {
          tanh(Layers.get(l)); //apply activation function
        }
      }
    } else
    {
      //output cannot be computed
    }

    return Layers.get(Layers.size()-1);
  }

  void RandomizeWeights(double weightsRange)
  {
    for (double[][] weight : Weights)
      for (int i = 0; i < weight.length; i++)
        for (int j = 0; j < weight[0].length; j++)
          weight[i][j] = ((Math.random()*2)-1) * weightsRange;
  }

  double MaxWeight()
  {
    double max = Double.MIN_VALUE;
    for (double[][] weight : Weights)
      for (int i = 0; i < weight.length; i++)
        for (int j = 0; j < weight[0].length; j++)
          if (weight[i][j] > max)
            max = weight[i][j];
    return max;
  }

  double MinWeight()
  {
    double min = Double.MAX_VALUE;
    for (double[][] weight : Weights)
      for (int i = 0; i < weight.length; i++)
        for (int j = 0; j < weight[0].length; j++)
          if (weight[i][j] < min)
            min = weight[i][j];
    return min;
  }

  void RandomizeBiases(double biases_range)
  {
    for (double[] Bias : Biases)
      for (int i = 0; i < Bias.length; i++)
        Bias[i] = ((Math.random()*2)-1) * biases_range;
  }

  double MaxBias()
  {
    double max = Double.MIN_VALUE;
    for (double[] Bias : Biases)
      for (int i = 0; i < Bias.length; i++)
        if (Bias[i] > max)
          Bias[i] = max;
    return max;
  }

  double MinBias()
  {
    double min = Double.MAX_VALUE;
    for (double[] Bias : Biases)
      for (int i = 0; i < Bias.length; i++)
        if (Bias[i] < min)
          Bias[i] = min;
    return min;
  }

  void Randomization( double weightsRange, double biasesRange)
  {
    RandomizeWeights(weightsRange);
    RandomizeBiases(biasesRange);
  }

  private double[] Mult(double[][] Matrix, double[] Vector)
  {
    double[] result = new double[Matrix.length]; //first dim of Matrix
    if (Matrix[0].length == Vector.length) //if matrix can be multiplied by vector
    {
      double temp;
      for (int i = 0; i<Matrix.length; i++)
      {
        temp = 0;
        for (int j = 0; j<Matrix[0].length; j++)
        {
          temp += Matrix[i][j] * Vector[j];
        }
        result[i] = temp;
      }
    } else
    {
      //result cannot be computed
    }
    return result;
  }

  private void Sum(double[] VectorA, double[] VectorB)
  {
    if (VectorA.length == VectorB.length)
    {
      for (int i=0; i<VectorA.length; i++ )
      {
        VectorA[i] = VectorA[i] + VectorB[i];
      }
    } else
    {
      //result cannot be computed
    }
  }

  /*ReLU activation*/
  private void ReLU(double[] arr)
  {
    for (int i = 0; i < arr.length; i++)
      if (arr[i]<=0) arr[i]=0;
  }
  private double ReLUDerivative(double x)
  {
    return x<0 ? 0 : 1;
  }

  /* tanh activation*/
  private void tanh(double[] arr)
  {
    for (int i = 0; i < arr.length; i++)
      arr[i] = Math.tanh(arr[i]);
  }
  private double tanhDerivative(double x)
  {
    double temp = Math.tanh(x);
    return 1-temp*temp;
  }

  public void showAt(float x, float y, float Width, float Height)
  {
    push();
    translate(x,y);
    show(Width, Height);
    pop();
  }

  public void show(float Width, float Height)
  {
    int horizontalDim = LayersDimensions.size();
    int verticalDim=0;
    float nodeSpanHorizontal;
    float nodeSpanVertical;
    float nodeDiagonal;
    float strokeRatio = 0.1;
    float gapRatio = 1; // gapWidth/nodeDiagonal
    ArrayList<PVector[]> Points = new ArrayList<PVector[]>();
    for (int i=0; i<LayersDimensions.size(); i++)
    {
      if (verticalDim < LayersDimensions.get(i))
        verticalDim = LayersDimensions.get(i);
    }
    float tempW = Width/(horizontalDim*(gapRatio+1)-gapRatio);
    float tempH = Height/(verticalDim*(gapRatio+1)-gapRatio);
    nodeDiagonal = tempW<tempH ? tempW : tempH;
    nodeSpanHorizontal = horizontalDim > 1 ? (Width - nodeDiagonal*horizontalDim)/(horizontalDim-1) : 0;
    nodeSpanVertical = verticalDim > 1 ? (Height - nodeDiagonal*verticalDim)/(verticalDim-1) : 0;
    for (int i = 0; i<LayersDimensions.size(); i++)
    {
      float horizontalBias = ((nodeDiagonal+nodeSpanHorizontal)*i)+nodeDiagonal/2;
      float verticalBias = (verticalDim-LayersDimensions.get(i))*(nodeDiagonal+nodeSpanVertical)/2 + nodeDiagonal/2;
      Points.add(new PVector[LayersDimensions.get(i)]);
      for (int j = 0; j<LayersDimensions.get(i); j++)
      {
        Points.get(i)[j] = new PVector(horizontalBias, verticalBias+((nodeDiagonal+nodeSpanVertical)*j));
      }
    }

    //draw

    //lines
    stroke(0);
    strokeWeight(nodeDiagonal*strokeRatio);
    double maxW = MaxWeight();
    double minW = MinWeight();
    for (int i = 0; i<LayersDimensions.size()-1; i++)
    {
      for (int j = 0; j<LayersDimensions.get(i); j++)
      {
        for (int k = 0; k<LayersDimensions.get(i+1); k++) 
        {
          float colour = (float)(((Weights.get(i)[j][k] - minW)*512)/(maxW-minW));
          stroke(512-colour, colour, colour<256 ?colour:(512-colour));
          line(Points.get(i)[j].x, Points.get(i)[j].y, Points.get(i+1)[k].x, Points.get(i+1)[k].y);
        }
      }
    }

    //circles
    strokeWeight(nodeDiagonal*strokeRatio);
    stroke(0);
    for (int i = 0; i<LayersDimensions.size(); i++)
    {
      for (int j = 0; j<LayersDimensions.get(i); j++)
      {
        circle(Points.get(i)[j].x, Points.get(i)[j].y, nodeDiagonal);
      }
    }
  }
}
