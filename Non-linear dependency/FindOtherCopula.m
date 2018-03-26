function [Copula] = FindOtherCopula(Y,NN,type)

%%%  Type = 'Gaussian' for Gaussian copula
%= 't' for t-copula
%= 'Clayton' for Clayton copula
%= 'Frank' for Frank copula
%= 'Independent' for independent copula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=[];
parfor i=1:size(Y,2)
    mc=ksdensity(Y(:,i),Y(:,i),'function','cdf');
    mc(find(mc<0.00001))=0.00001;
    mc(find(mc>0.99999))=0.99999;
    U=[U mc];
end

y=linspace(1/NN,1,NN);
Column=size(Y,2);

if Column==2
    i = 1:length(y); % indices into the vector x
    [i1,i2] = ndgrid(i); % all possible combinations of indices into x
    V = y([i2(:) i1(:)]); % 
elseif Column==3
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3] = ndgrid(i); % all possible combinations of indices into x
    V = y([i3(:) i2(:) i1(:)]); %
elseif Column==4
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4] = ndgrid(i); % all possible combinations of indices into x
    V = y([i4(:) i3(:) i2(:) i1(:)]); %
elseif Column==5
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5] = ndgrid(i); % all possible combinations of indices into x
    V = y([i5(:) i4(:) i3(:) i2(:) i1(:)]);
elseif Column==6
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6] = ndgrid(i); % all possible combinations of indices into x
    V = y([i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); % 
elseif Column==7
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6,i7] = ndgrid(i); % all possible combinations of indices into x
    V = y([i7(:) i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); % 
elseif Column==8
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6,i7,i8] = ndgrid(i); % all possible combinations of indices into x
    V = y([i8(:) i7(:) i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); %
elseif Column==9
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6,i7,i8,i9] = ndgrid(i); % all possible combinations of indices into x
    V = y([i9(:) i8(:) i7(:) i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); %
elseif Column==10
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10] = ndgrid(i); % all possible combinations of indices into x
    V = y([i10(:) i9(:) i8(:) i7(:) i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); %
elseif Column==11
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11] = ndgrid(i); % all possible combinations of indices into x
    V = y([i11(:) i10(:) i9(:) i8(:) i7(:) i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); %
elseif Column==12
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12] = ndgrid(i); % all possible combinations of indices into x
    V = y([i12(:) i11(:) i10(:) i9(:) i8(:) i7(:) i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); %
elseif Column==13
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13] = ndgrid(i); % all possible combinations of indices into x
    V = y([i13(:) i12(:) i11(:) i10(:) i9(:) i8(:) i7(:) i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); %   
elseif Column==14
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14] = ndgrid(i); % all possible combinations of indices into x
    V = y([i14(:) i13(:) i12(:) i11(:) i10(:) i9(:) i8(:) i7(:) i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); %   
elseif Column==15
    i = 1:length(y); % indices into the vector x
    [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15] = ndgrid(i); % all possible combinations of indices into x
    V = y([i15(:) i14(:) i13(:) i12(:) i11(:) i10(:) i9(:) i8(:) i7(:) i6(:) i5(:) i4(:) i3(:) i2(:) i1(:)]); %   
end


if strcmp(type,'Gaussian')==1
    Rho=copulafit('Gaussian',U);
    Y = copulacdf(type,V,Rho);
    Copula=reshape(Y,length(y)*ones(1,Column));
elseif strcmp(type,'Clayton')==1
    Rho=copulafit('Clayton',U);
    Y = copulacdf(type,V,Rho);
    Copula=reshape(Y,length(y)*ones(1,Column));
elseif strcmp(type,'Frank')==1
    Rho = copulafit('Frank',U);
    Y = copulacdf('Frank',V, Rho);
    Copula=reshape(Y,length(y)*ones(1,Column));
end
