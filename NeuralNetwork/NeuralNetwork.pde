int input = 1;
int output = 2;
int maxlayers = 5;
int maxnodes = 5;
neuralnetwork[][] nn = new neuralnetwork[maxlayers][maxnodes];
void setup() {
  size(1200, 800);
  surface.setResizable(true);

  for (int i=0; i<maxlayers; i++)
  {
    for (int j=0; j<maxnodes; j++)
    {
      nn[i][j] = new neuralnetwork(input, i+1, j+1, output);
      nn[i][j].RandomizeWeights(10);
    }
  }
}

void draw() {
  background(255);
  float nnWidth = (width*0.75)/maxlayers;
  float nnHeight = (height*0.75)/maxnodes;
  translate(nnWidth*0.15, nnHeight*0.15);
  for (int i=0; i<maxlayers; i++)
  {
    for (int j=0; j<maxnodes; j++)
    {
      nn[i][j].showAt((width/maxlayers)*i, (height/maxnodes)*j, nnWidth, nnHeight);
    }
  }
}

void keyPressed() {
  for (int i=0; i<maxlayers; i++)
  {
    for (int j=0; j<maxnodes; j++)
    {
      nn[i][j].RandomizeWeights(10);
    }
  }
}
