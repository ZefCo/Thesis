%  ANN to differentate Introns and Exons

%  Read Data and re-organize 
NumImages  = 3000;
ExonPath   = '../Images/Images_Q3_3K/Exon/';
IntronPath = '../Images/Images_Q3_3K/Intron/';

for nimag = 1:1:NumImages
    dtExon = imread (strcat(ExonPath, "Image", sprintf('%05d',nimag), ".png"));
    tmp    = dtExon(:, :, 1);
    Input(:, nimag)  = tmp(:);
    Target(1, nimag) = 0;

    dtIntron = imread (strcat(IntronPath, "Image", sprintf('%05d',nimag), ".png"));
    tmp      = dtIntron(:, :, 1);
    Input(:, NumImages+nimag)  = tmp(:);
    Target(1, NumImages+nimag) = 1;
end

%  Training
trainFcn             = 'trainscg';
hiddenLayerSize      = 10;
MolPatNet            = patternnet(hiddenLayerSize, trainFcn);
MolPatNet.divideFcn  = 'dividerand';  
MolPatNet.divideMode = 'sample';  

MolPatNet.divideParam.trainRatio = 50/100;
MolPatNet.divideParam.valRatio   = 15/100;
MolPatNet.divideParam.testRatio  = 35/100;
MolPatNet.performFcn             = 'crossentropy';

MolPatNet.plotFcns = {'plotperform', 'plottrainstate', 'ploterrhist', 'plotconfusion', 'plotroc'};
[MolPatNet,Trs]    = train(MolPatNet, double(Input), Target);
SimData            = MolPatNet (double (Input));
Diff               = gsubtract(Target, SimData);
performance        = perform(MolPatNet, Target, SimData);

trOut  = SimData(:, Trs.trainInd);
vOut   = SimData (:, Trs.valInd);
tsOut  = SimData (:, Trs.testInd);
trTarg = Target(:, Trs.trainInd);
vTarg  = Target (:, Trs.valInd);
tsTarg = Target (:,Trs.testInd);
figure, plotconfusion(trTarg, trOut, 'Train', vTarg, vOut, 'Validation', tsTarg, tsOut, 'Testing', Target,SimData,'All')
figure, plotroc(trTarg, trOut, 'Train', vTarg, vOut, 'Validation', tsTarg, tsOut, 'Testing', Target,SimData,'All')