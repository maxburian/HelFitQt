#include <qcustomplot.h>
#include "helfitqt.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	HelFitQt w;
	w.show();

	w.init();
	w.initialize_graphs();

	return a.exec();
}
