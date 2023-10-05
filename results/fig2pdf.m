function fig2pdf(fig,pdfname)
    if (~strcmp(pdfname(end-3:end),'.pdf'))
        tmpname = split(pdfname,'.');
        pdfname = strcat([tmpname(1),'.pdf']);
    end
    pdfname = strcat(['C:/Users/rydah/Documents/iwomp_poisson_23/figures/',pdfname]);
    exportgraphics(fig,pdfname,'BackgroundColor','none','Resolution',300)
    disp(['Printed figure to ',pdfname])
end