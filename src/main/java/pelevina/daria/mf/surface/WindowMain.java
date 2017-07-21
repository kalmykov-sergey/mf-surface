package pelevina.daria.mf.surface;

import javafx.animation.AnimationTimer;
import javafx.application.Application;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import javafx.stage.Stage;

import java.util.EnumSet;
import java.util.concurrent.*;

import static pelevina.daria.mf.surface.Constants.R;

public class WindowMain extends Application {

    private final Data data;
    private final Calculation calculation;
    private Canvas canvas;
    private GraphicsContext gc;
    private double padding = 100;
    private double scale = 1300;
    private CrossingFlag flag;
    protected AnimationTimer at = new AnimationTimer() {
        @Override
        public void handle(long now) {
            flag = calculation.next();
            drawContainer();
            drawSurface();
            if (!CrossingFlag.CONTINUE.contains(flag)) {
                at.stop();
            }
        }
    };

    public static void main(String[] args) {
        Application.launch(args);
    }

    public WindowMain() {
        data = new Data(4 * Constants.N);
        calculation = new Calculation(data);
    }

    @Override
    public void start(Stage primaryStage) throws Exception {
        primaryStage.setTitle("Расчет поверхности МЖ");
        Group root = new Group();
        double width = 2 * padding + getScaledLength(Constants.CONTAINER_RADIUS);
        double height = 2 * padding + getScaledLength(Constants.CONTAINER_HEIGHT);
        canvas = new Canvas(width, height);
        gc = canvas.getGraphicsContext2D();
        root.getChildren().add(canvas);
        Scene scene = new Scene(root);
        primaryStage.setScene(scene);
        scene.setOnKeyPressed(event -> {
            switch (event.getCode()) {
                case SPACE:
                    at.start();
                    break;
            }
        });
        primaryStage.show();
        flag = calculation.perform();
        drawContainer();
        drawSurface();


    }

    private double getScaledLength(double length) {
        return scale * length;
    }

    private double getX(double x) {
        return padding + scale * x;
    }

    private double getY(double y) {
        return canvas.getHeight() - padding - scale * y - getScaledLength(R);
    }

    public void drawContainer() {
        gc.setFill(Color.WHITE);
        gc.fillRect(0, 0, canvas.getWidth(), canvas.getHeight());
        gc.setStroke(Color.BLUE);
        gc.strokeLine(getX(0), getY(-R), getX(0), getY(Constants.CONTAINER_HEIGHT));
        gc.setStroke(Color.BLACK);
        gc.strokeOval(getX(-R), getY(R), getScaledLength(2 * R), getScaledLength(2 * R));
        gc.strokeLine(getX(0), getY(-R), getX(Constants.CONTAINER_RADIUS), getY(-R));
        gc.strokeLine(getX(Constants.CONTAINER_RADIUS), getY(-R), getX(Constants.CONTAINER_RADIUS), getY(Constants.CONTAINER_HEIGHT));
    }

    public void drawSurface() {
        gc.setStroke(Color.DARKRED);
        gc.setFill(Color.BLACK);
//        gc.setLineWidth(2);
        for (int j = 1; j <= data.lastIndex; j++) {
            Data.Point current = data.points[j];
            Data.Point previous = data.points[j - 1];
            gc.strokeLine(getX(previous.getX()), getY(previous.getY()), getX(current.getX()), getY(current.getY()));
            gc.fillOval(getX(previous.getX()) - 1, getY(previous.getY()) - 1, 2, 2);
            if (j == data.lastIndex) {
                gc.fillOval(getX(current.getX()) - 1, getY(current.getY()) - 1, 2, 2);
            }
        }
        gc.setLineWidth(1);
    }
}

