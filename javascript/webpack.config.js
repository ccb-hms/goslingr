var path = require('path');
var webpack = require('webpack');

module.exports = {
    mode: 'production',
    entry: './src/goslingr.js',
    output: {
        path: path.join(__dirname, '../inst/htmlwidgets'),
        filename: 'goslingr.js'
    },
    module: {
        rules: [
            {
                test: /\.jsx?$/,
                loader: 'babel-loader',
                options: {
                    presets: ['@babel/preset-env', '@babel/preset-react']
                }
            },
            // For CSS so that import "path/style.css"; works
            {
                test: /\.css$/,
                use: ['style-loader', 'css-loader']
            },
            { 
                test: /\.mjs$/, 
                include: /node_modules/, 
                type: 'javascript/auto' 
            }
        ]
    },
    stats: {
        colors: true
    },
    plugins: [
        new webpack.optimize.LimitChunkCountPlugin({
            maxChunks: 1
        })
    ]
};