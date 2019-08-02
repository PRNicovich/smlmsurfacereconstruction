function gauss = gauss_2D(parameters, domain)

Amp = parameters(1);
center_x = parameters(2);
center_y = parameters(3);
sigma_x = parameters(4);
sigma_y = parameters(5);
Background = parameters(6);

[D_x, D_y] = meshgrid(domain);

gauss = Amp*exp(-((((D_x - center_x).^2)./(2*sigma_x.^2)) + (((D_y - center_y).^2)./(2*sigma_y.^2)))) + Background;