using GrpcServer.Interceptors;
using GrpcServer.Services;

var builder = WebApplication.CreateBuilder(args);

// Add services to the container.
var serviceProvider = builder.Services.BuildServiceProvider();
var logger = serviceProvider.GetService<ILogger<LoggingInterceptor>>();
builder.Services.AddGrpc(options =>
{
    options.Interceptors.Add(typeof(LoggingInterceptor), logger);
});

var app = builder.Build();

// Configure the HTTP request pipeline.
app.MapGrpcService<GreeterService>();
app.MapGet("/", () => "Communication with gRPC endpoints must be made through a gRPC client. To learn how to create a client, visit: https://go.microsoft.com/fwlink/?linkid=2086909");

app.Run();
