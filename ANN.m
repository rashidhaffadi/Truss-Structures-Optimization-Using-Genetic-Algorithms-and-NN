function z=ANN(x)
    filename='./model/net';
    model=load(filename);
    z=model.net(x');
end

